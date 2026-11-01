import { randomUUID } from 'node:crypto'
import { execFileSync } from 'node:child_process'
import { mkdtemp, mkdir, readFile, rename, rm, writeFile } from 'node:fs/promises'
import { basename, dirname, join, relative, resolve } from 'node:path'
import { fileURLToPath, pathToFileURL } from 'node:url'
import { gzipSync } from 'node:zlib'
import { geneBucket } from '../src/data/geneBucket.mjs'
import { exists, files, hash, sourceInventory } from './release-sources.mjs'
import { exportAtlas } from './export.mjs'
import { exportReports, studyReportNames } from './export-reports.mjs'
import { releaseInventory } from './release-inventory.mjs'
import { compareAnalysisMetadata, releaseContext } from './release-provenance.mjs'
import { stampAssets } from './release-assets.mjs'
import { verifyRelease } from './verify-release.mjs'
import { pantherFilename, pantherRelease, pantherSha256 } from './panther-release.mjs'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
function revision() {
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
  const input = resolve(inputPath),
    target = resolve(targetPath),
    parent = dirname(target)
  await mkdir(parent, { recursive: true })
  const stage = await mkdtemp(join(parent, '.atlas-stage-'))
  const backup = join(parent, `.atlas-backup-${randomUUID()}`)
  let moved = false
  try {
    const before = await sourceInventory(input)
    const counts = await exportAtlas(input, join(stage, 'de-index.json'))
    const pathwayDirectory = resolve(options.pathwayDirectory || join(root, 'public/pathways'))
    const env = { ...process.env, SITE_DATA_DIR: stage, SITE_PATHWAY_DIR: pathwayDirectory }
    const pathwayInventory = async () => {
      const found = []
      for (const path of await files(pathwayDirectory))
        if (path.endsWith('.gmt') || path.endsWith('pathway-labels.json'))
          found.push({ path: relative(pathwayDirectory, path), sha256: hash(await readFile(path)) })
      return found
    }
    const pathwayFiles = await pathwayInventory()
    execFileSync(process.execPath, [join(root, 'scripts/build-pathway-index.mjs')], { env, stdio: 'inherit' })
    const panther =
      options.pantherPath || join(root, 'data-inputs/reference-data/panther', pantherRelease, pantherFilename)
    if ((await exists(panther)) && (await exists(join(stage, 'pipeline-trajectory-index.json')))) {
      if (hash(await readFile(panther)) !== pantherSha256)
        throw new Error('PANTHER source checksum does not match the pinned release')
      execFileSync(process.execPath, [join(root, 'scripts/build-gene-metadata-index.mjs')], {
        env,
        stdio: 'inherit',
      })
    }
    const trajectorySource = join(input, '..', 'combined', 'trajectory')
    const trajectoryCsv = await exists(join(trajectorySource, 'cross_cell_type_trajectory_clusters.csv'))
    const trajectoryMeta = await exists(join(trajectorySource, 'analysis_metadata.json'))
    if (trajectoryCsv !== trajectoryMeta)
      throw new Error('Trajectory inputs are incomplete: CSV and metadata must both be present')
    if (trajectoryCsv && counts.trajectory.input > 0 && counts.trajectory.exported === 0)
      throw new Error('Trajectory input produced no valid records; inspect age-bin support and input schema')
    const after = await sourceInventory(input)
    if (JSON.stringify(pathwayFiles) !== JSON.stringify(await pathwayInventory()))
      throw new Error('Pathway sources changed during export; retry from a stable snapshot')
    if (JSON.stringify(before) !== JSON.stringify(after))
      throw new Error('Analysis inputs changed during export; retry from a stable snapshot')
    const id = randomUUID()
    const trajectoryPath = join(stage, 'pipeline-trajectory-index.json')
    if (await exists(trajectoryPath)) {
      const trajectory = JSON.parse(await readFile(trajectoryPath, 'utf8'))
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
    for (const [cellIndex, cell] of catalogue.cellTypes.entries()) {
      for (const [modelIndex, model] of cell.models.entries()) {
        const asset = `de-models/${cellIndex}-${modelIndex}.json`
        await mkdir(join(stage, 'de-models'), { recursive: true })
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
              source: catalogue.source,
              cellTypes: [{ ...cell, models: [{ ...model, records }] }],
            }),
          )
          recordBuckets[bucket] = path
        }
        await writeFile(
          join(stage, asset),
          JSON.stringify({
            schemaVersion: 1,
            source: catalogue.source,
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
      }
    }
    await writeFile(dePath, JSON.stringify(catalogue))
    const assets = await stampAssets(stage, id)
    for (const asset of assets) {
      const source = await readFile(join(stage, asset.path))
      const compressed = gzipSync(source, { level: 9 })
      if (compressed.length + 128 >= source.length) continue
      asset.gzipPath = `${asset.path}.gz.bin`
      asset.gzipBytes = compressed.length
      asset.gzipSha256 = hash(compressed)
      await writeFile(join(stage, asset.gzipPath), compressed)
    }
    const de = JSON.parse(await readFile(join(stage, 'de-index.json'), 'utf8'))
    let reportSources = []
    if (options.includeReports !== false) {
      const snapshotReports = async () => {
        const names = [
          ...(await studyReportNames(input)),
          ...de.cellTypes.map((cell) => `${cell.id}/differential_expression/report.html`),
          '../combined/trajectory/report.html',
        ]
        const snapshot = []
        for (const name of names) {
          const path = resolve(input, name)
          if (await exists(path)) snapshot.push({ path: name, sha256: hash(await readFile(path)) })
        }
        return snapshot
      }
      reportSources = await snapshotReports()
      assets.push(...(await exportReports(input, stage, id, de.cellTypes)))
      if (JSON.stringify(reportSources) !== JSON.stringify(await snapshotReports()))
        throw new Error('HTML reports changed during export; retry from a stable snapshot')
    }
    const { inventory, analysis } = await releaseInventory(input, de)
    const provenance = await compareAnalysisMetadata(input, inventory)
    const studiesPath = resolve(options.studiesConfigPath || join(root, '..', 'config/studies.json'))
    const context = await releaseContext(input, de, studiesPath)
    const referencesAfter = await releaseContext(input, de, studiesPath)
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
    await verifyRelease(stage)
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
