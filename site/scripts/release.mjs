import { randomUUID } from 'node:crypto'
import { execFileSync } from 'node:child_process'
import { copyFile, mkdtemp, mkdir, readFile, rename, rm, writeFile } from 'node:fs/promises'
import { basename, dirname, join, relative, resolve } from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'
import { gzip, gzipSync } from 'node:zlib'
import { promisify } from 'node:util'
import { geneBucket } from '../src/data/geneBucket.mjs'
import { exists, files, hash, hashFile, sourceInventory } from './release-sources.mjs'
import { exportAtlas } from './export.mjs'
import { exportReports, studyReportNames } from './export-reports.mjs'
import { releaseInventory } from './release-inventory.mjs'
import { compareAnalysisMetadata, releaseContext } from './release-provenance.mjs'
import { stampAssets } from './release-assets.mjs'
import { verifyRelease } from './verify-release.mjs'
import { pantherFilename, pantherRelease, pantherSha256 } from './panther-release.mjs'
import { exportSingleCell } from './export-single-cell.mjs'
import { exportClocks } from './export-clocks.mjs'
import { addStudySupport } from './export-study-support.mjs'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const gzipAsync = promisify(gzip)
function differentialPreviewSvg(records, label) {
  const width = 900, height = 520, left = 76, right = 24, top = 48, bottom = 58
  const safeLabel = label.replaceAll('&', '&amp;').replaceAll('<', '&lt;').replaceAll('>', '&gt;').replaceAll('"', '&quot;').replaceAll("'", '&apos;')
  const finite = records.filter((record) => Number.isFinite(record.effect) && Number.isFinite(record.padj))
  if (!finite.length) return null
  const nonSignificant = finite.filter((record) => record.padj >= 0.05)
  const significant = finite.filter((record) => record.padj < 0.05)
  const sample = (rows, maximum) => rows.length <= maximum
    ? rows
    : rows.filter((_, index) => index % Math.ceil(rows.length / maximum) === 0).slice(0, maximum)
  const visible = [...sample(nonSignificant, 2400), ...sample(significant, 3600)]
  const magnitudes = finite.map((record) => Math.abs(record.effect)).sort((a, b) => a - b)
  const xLimit = Math.max(0.5, magnitudes[Math.floor((magnitudes.length - 1) * 0.98)])
  const maxScore = finite.reduce((maximum, record) => Math.max(maximum, -Math.log10(Math.max(record.padj, 1e-15))), 0)
  const yLimit = Math.max(2, Math.min(15, maxScore))
  const x = (value) => left + ((Math.max(-xLimit, Math.min(xLimit, value)) + xLimit) / (2 * xLimit)) * (width - left - right)
  const y = (value) => height - bottom - (Math.min(yLimit, value) / yLimit) * (height - top - bottom)
  const points = visible.map((record) => {
    const score = -Math.log10(Math.max(record.padj, 1e-15))
    const significant = record.padj < 0.05
    return `<circle cx="${x(record.effect).toFixed(1)}" cy="${y(score).toFixed(1)}" r="2.1" fill="${significant ? '#3b8c78' : '#9ca6a3'}" fill-opacity="${significant ? '0.72' : '0.32'}"/>`
  }).join('')
  const threshold = -Math.log10(0.05)
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-labelledby="title desc"><title id="title">Differential expression preview for ${safeLabel}</title><desc id="desc">Static preview of the combined age model volcano plot. Interactive results are loading.</desc><rect width="100%" height="100%" fill="#fbfaf7"/><text x="${left}" y="28" font-family="sans-serif" font-size="17" fill="#263633">${safeLabel} · Combined age model</text><path d="M${left} ${top}V${height-bottom}H${width-right}" fill="none" stroke="#6d7773" stroke-width="1"/><path d="M${left} ${y(threshold).toFixed(1)}H${width-right}" stroke="#9d7d45" stroke-dasharray="5 5" stroke-width="1"/>${points}<text x="${(left + width - right) / 2}" y="${height-16}" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#44524e">Effect (log₂ fold change)</text><text transform="translate(18 ${(top + height - bottom) / 2}) rotate(-90)" text-anchor="middle" font-family="sans-serif" font-size="13" fill="#44524e">−log₁₀ adjusted p-value</text><text x="${width-right-8}" y="${(y(threshold)-6).toFixed(1)}" text-anchor="end" font-family="sans-serif" font-size="11" fill="#806538">adjusted p = 0.05</text></svg>\n`
}
const formatElapsed = (milliseconds) => {
  const seconds = Math.floor(milliseconds / 1000)
  return `${String(Math.floor(seconds / 60)).padStart(2, '0')}:${String(seconds % 60).padStart(2, '0')}`
}
const formatBytes = (bytes) => `${(bytes / 1024 / 1024).toFixed(1)} MiB`
function revision() {
  const configuredRevision = process.env.SITE_EXPORTER_REVISION?.trim()
  if (configuredRevision) return configuredRevision
  try {
    return execFileSync('git', ['rev-parse', 'HEAD'], {
      cwd: root,
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'ignore'],
    }).trim()
  } catch {
    return null
  }
}
export async function createRelease(inputPath, targetPath, options = {}) {
  const startedAt = Date.now()
  const progress = (message) =>
    console.log(`[site-export +${formatElapsed(Date.now() - startedAt)}] ${message}`)
  const timed = async (label, action) => {
    const phaseStartedAt = Date.now()
    progress(`Starting ${label}`)
    try {
      const result = await action()
      progress(`Finished ${label} in ${formatElapsed(Date.now() - phaseStartedAt)}`)
      return result
    } catch (error) {
      progress(`Failed ${label} after ${formatElapsed(Date.now() - phaseStartedAt)}`)
      throw error
    }
  }
  const slowHash = (name, bytesRead) => progress(`Still hashing ${name} (${formatBytes(bytesRead)} read)`)
  const input = resolve(inputPath),
    target = resolve(targetPath),
    parent = dirname(target)
  await mkdir(parent, { recursive: true })
  const stage = await mkdtemp(join(parent, '.atlas-stage-'))
  const backup = join(parent, `.atlas-backup-${randomUUID()}`)
  let moved = false
  try {
    progress(`Preparing release from ${input}`)
    const before = await timed('checking and hashing analysis inputs', () => sourceInventory(input, slowHash))
    const counts = await timed('exporting differential-expression and summary data', () => exportAtlas(input, join(stage, 'de-index.json'), {
      writeModel: async (cell, model, cellIndex, modelIndex, sourceFile) => {
        const asset = `de-models/${cellIndex}-${modelIndex}.json`
        await mkdir(join(stage, 'de-models'), { recursive: true })
        if (model.model !== 'per_study' && model.model !== 'developmental') {
          const downloadAsset = `de-downloads/${cell.id}/${model.model}-${model.covariate}.csv`
          const downloadPath = join(stage, downloadAsset)
          await mkdir(dirname(downloadPath), { recursive: true })
          await copyFile(sourceFile, downloadPath)
          model.downloadAsset = downloadAsset
        }
        if (model.model === 'merged' && model.covariate === 'age') {
          const preview = differentialPreviewSvg(model.records, cell.label)
          if (preview) {
            const previewAsset = `de-previews/${cell.id}.svg`
            const previewPath = join(stage, previewAsset)
            await mkdir(dirname(previewPath), { recursive: true })
            await writeFile(previewPath, preview)
            cell.previewAsset = previewAsset
          }
        }
        const allRecords = model.records
        const otherBuckets = new Map()
        for (const record of allRecords) {
          if (record.padj !== null && record.padj < 0.05) continue
          const bucket = geneBucket(record.gene)
          if (!otherBuckets.has(bucket)) otherBuckets.set(bucket, [])
          otherBuckets.get(bucket).push(record)
        }
        const recordBuckets = {}
        for (const [bucket, records] of otherBuckets) {
          const path = `de-models/${cellIndex}-${modelIndex}-${bucket}.json`
          await writeFile(
            join(stage, path),
            JSON.stringify({
              schemaVersion: 1,
              source: 'pseudobulk differential expression',
              cellTypes: [{ ...cell, models: [{ ...model, records }] }],
            }),
          )
          recordBuckets[bucket] = path
        }
        await writeFile(
          join(stage, asset),
          JSON.stringify({
            schemaVersion: 1,
            source: 'pseudobulk differential expression',
            cellTypes: [
              {
                ...cell,
                models: [
                  {
                    ...model,
                    records: allRecords.filter((record) => record.padj !== null && record.padj < 0.05),
                  },
                ],
              },
            ],
          }),
        )
        model.records = []
        model.recordsAsset = asset
        model.recordBuckets = recordBuckets
      },
    }))
    const pathwayDirectory = resolve(options.pathwayDirectory || join(root, 'public/pathways'))
    const env = { ...process.env, SITE_DATA_DIR: stage, SITE_PATHWAY_DIR: pathwayDirectory }
    const pathwayInventory = async () => {
      const found = []
      for (const path of await files(pathwayDirectory))
        if (path.endsWith('.gmt') || path.endsWith('pathway-labels.json'))
          found.push({ path: relative(pathwayDirectory, path), sha256: await hashFile(path) })
      return found
    }
    const pathwayFiles = await timed('checking pathway source files', pathwayInventory)
    await timed('building the pathway index', async () =>
      execFileSync(process.execPath, [join(root, 'scripts/build-pathway-index.mjs')], { env, stdio: 'inherit' }),
    )
    const panther =
      options.pantherPath || join(root, 'data-inputs/reference-data/panther', pantherRelease, pantherFilename)
    if ((await exists(panther)) && (await exists(join(stage, 'pipeline-trajectory-index.json')))) {
      if ((await hashFile(panther)) !== pantherSha256)
        throw new Error('PANTHER source checksum does not match the pinned release')
      await timed('building the gene metadata index', async () =>
        execFileSync(process.execPath, [join(root, 'scripts/build-gene-metadata-index.mjs')], {
          env: { ...env, SITE_PANTHER_PATH: resolve(panther) },
          stdio: 'inherit',
        }),
      )
    }
    const trajectorySource = join(input, '..', 'combined', 'trajectory')
    const trajectoryCsv = await exists(join(trajectorySource, 'cross_cell_type_trajectory_clusters.csv'))
    const trajectoryMeta = await exists(join(trajectorySource, 'analysis_metadata.json'))
    if (trajectoryCsv !== trajectoryMeta)
      throw new Error('Trajectory inputs are incomplete: CSV and metadata must both be present')
    if (trajectoryCsv && counts.trajectory.input > 0 && counts.trajectory.exported === 0)
      throw new Error('Trajectory input produced no valid records; inspect age-bin support and input schema')
    progress('Rechecking analysis inputs for changes during export')
    const after = await sourceInventory(input, slowHash)
    if (JSON.stringify(pathwayFiles) !== JSON.stringify(await pathwayInventory()))
      throw new Error('Pathway sources changed during export; retry from a stable snapshot')
    if (JSON.stringify(before) !== JSON.stringify(after))
      throw new Error('Analysis inputs changed during export; retry from a stable snapshot')
    const id = randomUUID()
    const trajectoryPath = join(stage, 'pipeline-trajectory-index.json')
    if (await exists(trajectoryPath)) {
      const trajectory = JSON.parse(await readFile(trajectoryPath, 'utf8'))
      const mapBuckets = new Map()
      for (const record of trajectory.records) {
        const bucket = geneBucket(record.gene)
        if (!mapBuckets.has(bucket)) mapBuckets.set(bucket, [])
        mapBuckets.get(bucket).push(record)
      }
      const trajectoryMap = {
        schemaVersion: 1,
        release: trajectory.release,
        warning: trajectory.warning,
        description: trajectory.description,
        records: trajectory.records
          .filter((record) => record.adjustedPValue < 0.05)
          .map(({ id, x, y, cellType, tissue, gene, adjustedPValue, cluster }) =>
            ({ id, x, y, cellType, tissue, gene, adjustedPValue, cluster })),
        tissues: trajectory.tissues,
        cellTypes: trajectory.cellTypes,
        genes: trajectory.genes,
        geneBuckets: {},
      }
      if (mapBuckets.size) await mkdir(join(stage, 'trajectory-details'), { recursive: true })
      for (const [bucket, records] of mapBuckets) {
        const path = `trajectory-details/${bucket}.json`
        await writeFile(
          join(stage, path),
          JSON.stringify({
            schemaVersion: 1,
            release: trajectory.release,
            records,
            genes: [...new Set(records.map((record) => record.gene))],
            cellTypes: trajectory.cellTypes,
            tissues: trajectory.tissues,
          }),
        )
        trajectoryMap.geneBuckets[bucket] = path
      }
      await writeFile(join(stage, 'pipeline-trajectory-map.json'), JSON.stringify(trajectoryMap))
      const buckets = new Map()
      for (const record of trajectory.records) {
        if (record.adjustedPValue < 0.05) continue
        const bucket = geneBucket(record.gene)
        if (!buckets.has(bucket)) buckets.set(bucket, [])
        buckets.get(bucket).push(record)
      }
      trajectory.geneBuckets = {}
      if (buckets.size) await mkdir(join(stage, 'trajectory-extras'), { recursive: true })
      for (const [bucket, records] of buckets) {
        const path = `trajectory-extras/${bucket}.json`
        await writeFile(
          join(stage, path),
          JSON.stringify({
            schemaVersion: 1,
            release: trajectory.release,
            records,
            genes: [...new Set(records.map((record) => record.gene))],
            cellTypes: trajectory.cellTypes,
            tissues: trajectory.tissues,
          }),
        )
        trajectory.geneBuckets[bucket] = path
      }
      trajectory.records = trajectory.records.filter((record) => record.adjustedPValue < 0.05)
      trajectory.description +=
        ' The default map shows adjusted p < 0.05; requested genes may include other retained trajectories.'
      await writeFile(trajectoryPath, JSON.stringify(trajectory))
    }
    // Keep the model catalogue small; fetch gene records only for the selected model.
    const dePath = join(stage, 'de-index.json')
    const catalogue = JSON.parse(await readFile(dePath, 'utf8'))
    for (const cell of catalogue.cellTypes) {
      if (!cell.trajectoryAsset) continue
      const previewPath = join(stage, cell.trajectoryAsset)
      const preview = JSON.parse(await readFile(previewPath, 'utf8'))
      const buckets = new Map()
      for (const [gene, bins] of Object.entries(preview.genes)) {
        const bucket = geneBucket(gene)
        if (!buckets.has(bucket)) buckets.set(bucket, {})
        buckets.get(bucket)[gene] = bins
      }
      const prefix = cell.trajectoryAsset.slice(0, -5)
      await mkdir(join(stage, prefix), { recursive: true })
      preview.geneBuckets = {}
      for (const [bucket, genes] of buckets) {
        const path = `${prefix}/${bucket}.json`
        await writeFile(
          join(stage, path),
          JSON.stringify({
            schemaVersion: 1,
            referenceBin: preview.referenceBin,
            design: preview.design,
            genes,
          }),
        )
        preview.geneBuckets[bucket] = path
      }
      preview.genes = {}
      await writeFile(previewPath, JSON.stringify(preview))
    }
    await writeFile(dePath, JSON.stringify(catalogue))
    await timed('preparing study support counts', () => addStudySupport(stage))
    await timed('preparing single-cell summaries', () => exportSingleCell(input, stage))
    await timed('preparing transcriptomic clock results', () => exportClocks(input, stage, catalogue))
    const assets = await timed('validating and stamping generated assets', () =>
      stampAssets(stage, id, (processed, total, valid) =>
        progress(`Validated ${processed}/${total} generated files (${valid} published assets)`),
      ),
    )
    let compressedCount = 0
    for (const [index, asset] of assets.entries()) {
      const source = await readFile(join(stage, asset.path))
      if (source.length >= 16 * 1024 * 1024)
        progress(`Compressing large asset ${index + 1}/${assets.length}: ${asset.path} (${formatBytes(source.length)})`)
      const compression = gzipAsync(source, { level: 9 })
      let heartbeat
      if (source.length >= 16 * 1024 * 1024)
        heartbeat = setTimeout(
          () => progress(`Still compressing ${asset.path} (${formatBytes(source.length)})`),
          30_000,
        )
      let compressed
      try {
        compressed = await compression
      } finally {
        if (heartbeat) clearTimeout(heartbeat)
      }
      if (compressed.length + 128 >= source.length) continue
      asset.gzipPath = `${asset.path}.gz.bin`
      asset.gzipBytes = compressed.length
      asset.gzipSha256 = hash(compressed)
      await writeFile(join(stage, asset.gzipPath), compressed)
      compressedCount += 1
      if (compressedCount % 25 === 0)
        progress(`Compressed ${compressedCount} assets so far`)
    }
    progress(`Compression complete: ${compressedCount} assets received gzip copies`)
    const de = JSON.parse(await readFile(join(stage, 'de-index.json'), 'utf8'))
    let reportSources = []
    if (options.includeReports !== false) {
      const snapshotReports = async () => {
        const names = [
          '../../run_manifest.json',
          ...(await studyReportNames(input)),
          ...(await studyReportNames(input)).map((path) => path.replace(/report\.html$/, 'executed.ipynb')),
          ...de.cellTypes.map((cell) => `${cell.id}/differential_expression/report.html`),
          '../combined/trajectory/report.html',
          'age_clock_cross_type/report.html',
        ]
        const snapshot = []
        for (const name of names) {
          const path = resolve(input, name)
          if (await exists(path))
            snapshot.push({
              path: name,
              sha256: await hashFile(path, (bytesRead) =>
                progress(`Still hashing report source ${name} (${formatBytes(bytesRead)} read)`),
              ),
            })
        }
        return snapshot
      }
      reportSources = await timed('checking report source files', snapshotReports)
      assets.push(
        ...(await timed('copying and indexing reports', () =>
          exportReports(input, stage, id, de.cellTypes, progress),
        )),
      )
      if (JSON.stringify(reportSources) !== JSON.stringify(await snapshotReports()))
        throw new Error('HTML reports changed during export; retry from a stable snapshot')
    }
    const { inventory, analysis } = await timed('building the analysis inventory', () => releaseInventory(input, de))
    const provenance = await timed('checking analysis provenance', () => compareAnalysisMetadata(input, inventory))
    const studiesPath = resolve(options.studiesConfigPath || join(root, '..', 'config/studies.json'))
    const context = await timed('collecting release context', () => releaseContext(input, de, studiesPath, assets))
    const referencesAfter = await releaseContext(input, de, studiesPath, assets)
    if (
      JSON.stringify(context) !== JSON.stringify(referencesAfter) ||
      JSON.stringify(before) !== JSON.stringify(await sourceInventory(input))
    )
      throw new Error('Analysis or reference metadata changed during export; retry from a stable snapshot')
    const manifest = {
      schemaVersion: 1,
      releaseId: id,
      generatedAt: new Date().toISOString(),
      exporterRevision: revision(),
      source: {
        input: basename(input),
        files: before,
        note: 'File hashes identify observed inputs; they do not prove that the inputs originate from one analysis run. Exporter revision is not necessarily the analysis revision.',
      },
      analysis,
      provenance,
      context,
      assets,
      inventory,
      counts,
      pathwaySources: pathwayFiles,
      reportSources,
      availability: {
        trajectory: assets.some((a) => a.kind === 'trajectory'),
        composition: assets.some((a) => a.kind === 'composition'),
        geneMetadata: assets.some((a) => a.kind === 'metadata'),
        clocks: assets.some((a) => a.kind === 'clock'),
        pathways: assets.some((a) => a.kind === 'pathways'),
        reports: assets.some((a) => a.kind === 'report'),
      },
      limitations: [
        'Pipeline run identity and source cohort licences must be verified separately.',
        'Trajectory records are the subset selected for clustering.',
      ],
    }
    for (const [name, count] of Object.entries(counts))
      if (count.input !== count.exported)
        console.warn(
          `${name}: ${count.input - count.exported} rows excluded by export validation (see release.json)`,
        )
    if (!analysis) console.warn('DE analysis manifest unavailable; release provenance is incomplete')
    const manifestBytes = Buffer.from(JSON.stringify(manifest, null, 2) + '\n')
    await writeFile(join(stage, 'release.json'), manifestBytes)
    await writeFile(join(stage, 'release.json.gz.bin'), gzipSync(manifestBytes, { level: 9 }))
    await timed('verifying the staged release', () => verifyRelease(stage))
    progress('Publishing the validated release')
    if (await exists(target)) {
      await rename(target, backup)
      moved = true
    }
    try {
      await rename(stage, target)
    } catch (e) {
      if (moved) {
        await rename(backup, target)
        moved = false
      }
      throw e
    }
    if (moved) await rm(backup, { recursive: true, force: true })
    progress(`Release ${id} published with ${assets.length} assets in ${formatElapsed(Date.now() - startedAt)}`)
    return manifest
  } finally {
    await rm(stage, { recursive: true, force: true })
  }
}
if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href)
  await createRelease(
    process.argv[2] || 'output/analysis/per_type',
    process.argv[3] || join(root, 'public/data'),
  )
