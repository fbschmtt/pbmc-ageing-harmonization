import { createHash } from 'node:crypto'
import { readFile, stat } from 'node:fs/promises'
import { isAbsolute, relative, resolve } from 'node:path'
import { pathToFileURL } from 'node:url'
import { gunzipSync } from 'node:zlib'
import { validAsset } from '../src/data/contracts.mjs'
import { clockModelAssetPath } from './clock-assets.mjs'
import { hashFile } from './release-sources.mjs'

export async function verifyRelease(directory) {
  const root = resolve(directory)
  const startedAt = Date.now()
  const progress = (message) => {
    const seconds = Math.floor((Date.now() - startedAt) / 1000)
    const elapsed = `${String(Math.floor(seconds / 60)).padStart(2, '0')}:${String(seconds % 60).padStart(2, '0')}`
    console.log(`[site-verify +${elapsed}] ${message}`)
  }
  progress(`Reading release manifest from ${root}`)
  const manifestBytes = await readFile(resolve(root, 'release.json'))
  const manifest = JSON.parse(manifestBytes.toString('utf8'))
  try {
    const compressed = await readFile(resolve(root, 'release.json.gz.bin'))
    if (!gunzipSync(compressed).equals(manifestBytes)) throw new Error('Compressed manifest mismatch')
  } catch (error) {
    if (error.code !== 'ENOENT') throw error
  }
  if (
    manifest.schemaVersion !== 1 ||
    typeof manifest.releaseId !== 'string' ||
    !Array.isArray(manifest.assets)
  )
    throw new Error('Invalid release manifest')
  const seen = new Set()
  const parsed = new Map()
  let lastProgress = Date.now()
  progress(`Checking ${manifest.assets.length} assets`)
  for (const [index, asset] of manifest.assets.entries()) {
    if (
      index === 0 ||
      index % 100 === 0 ||
      Date.now() - lastProgress >= 15_000 ||
      asset.bytes >= 50 * 1024 * 1024
    ) {
      progress(`Checking asset ${index + 1}/${manifest.assets.length}: ${asset.path}`)
      lastProgress = Date.now()
    }
    if (typeof asset.path !== 'string' || seen.has(asset.path))
      throw new Error('Invalid or duplicate asset path')
    seen.add(asset.path)
    const path = resolve(root, asset.path),
      rel = relative(root, path)
    if (rel.startsWith('..') || isAbsolute(rel)) throw new Error('Unsafe asset path')
    const streamReport =
      (asset.kind === 'report' || asset.kind === 'report_index') && asset.gzipPath === undefined
    const bytes = streamReport ? null : await readFile(path)
    if (
      (streamReport ? (await stat(path)).size : bytes.length) !== asset.bytes ||
      (streamReport
        ? await hashFile(path, (bytesRead) =>
            progress(`Still hashing ${asset.path} (${(bytesRead / 1024 / 1024).toFixed(1)} MiB read)`),
          )
        : createHash('sha256').update(bytes).digest('hex')) !==
        asset.sha256
    )
      throw new Error(`Asset checksum mismatch: ${asset.path}`)
    if (asset.gzipPath !== undefined) {
      if (asset.gzipPath !== `${asset.path}.gz.bin`) throw new Error('Invalid compressed asset path')
      const compressed = await readFile(resolve(root, asset.gzipPath))
      if (
        compressed.length !== asset.gzipBytes ||
        createHash('sha256').update(compressed).digest('hex') !== asset.gzipSha256 ||
        !gunzipSync(compressed).equals(bytes)
      )
        throw new Error(`Compressed asset mismatch: ${asset.path}`)
    }
    if (asset.kind === 'report' || asset.kind === 'report_index') {
      if (
        !asset.path.startsWith(`reports/${manifest.releaseId}/`) ||
        !(asset.kind === 'report'
          ? asset.path.endsWith('/report.html')
          : asset.path === `reports/${manifest.releaseId}/index.html`) ||
        typeof asset.label !== 'string'
      )
        throw new Error(`Invalid report asset: ${asset.path}`)
      continue
    }
    if (asset.kind === 'study_umap' || asset.kind === 'study_heatmap') {
      if (
        !asset.path.startsWith(`images/${manifest.releaseId}/studies/`) ||
        !asset.path.endsWith('.png') ||
        typeof asset.study !== 'string' ||
        !bytes.subarray(0, 8).equals(Buffer.from('89504e470d0a1a0a', 'hex'))
      )
        throw new Error(`Invalid study figure: ${asset.path}`)
      continue
    }
    if (asset.kind === 'de_download') {
      if (!asset.path.startsWith('de-downloads/') || !asset.path.endsWith('.csv'))
        throw new Error(`Invalid combined-model download: ${asset.path}`)
      continue
    }
    if (asset.kind === 'de_preview') {
      if (!asset.path.startsWith('de-previews/') || !asset.path.endsWith('.svg') ||
          !bytes.toString('utf8').startsWith('<svg xmlns="http://www.w3.org/2000/svg"'))
        throw new Error(`Invalid differential-expression preview: ${asset.path}`)
      continue
    }
    if (asset.kind === 'trajectory_preview') {
      if (asset.path !== 'pipeline-trajectory-preview.svg' ||
          !bytes.toString('utf8').startsWith('<svg xmlns="http://www.w3.org/2000/svg"'))
        throw new Error('Invalid trajectory map preview image')
      continue
    }
    const value = JSON.parse(bytes.toString('utf8'))
    if (value.releaseId !== manifest.releaseId || !validAsset(asset.kind, value))
      throw new Error(`Invalid release asset: ${asset.path}`)
    // Reference checks need only chunk identity and metadata, not their gene payloads.
    if (asset.path.startsWith('de-models/')) {
      for (const cell of value.cellTypes)
        for (const model of cell.models)
          model.records = model.records
            .filter((record) => record.padj === null || record.padj >= 0.05)
            .slice(0, 1)
    } else if (asset.path.startsWith('age-previews/')) {
      value.genes = {}
    } else if (asset.path.startsWith('trajectory-extras/')) {
      value.records = []
      value.genes = []
    }
    parsed.set(asset.path, value)
  }
  progress('Checking asset references and cell types')
  const singleCellIndex = parsed.get('single-cell-index.json')
  if (singleCellIndex) {
    const ids = new Set()
    for (const cell of singleCellIndex.cellTypes) {
      const details = parsed.get(cell.path)
      if (ids.has(cell.id) || cell.path !== `single-cell/${cell.id}.json` || !details || details.cellType !== cell.label)
        throw new Error(`Invalid single-cell reference: ${cell.id}`)
      ids.add(cell.id)
    }
  }
  if (!seen.has('de-index.json')) throw new Error('Missing DE asset')
  const de = parsed.get('de-index.json')
  const clocks = parsed.get('clock-data.json')
  if (!clocks && manifest.assets.some((asset) => asset.kind === 'clockModel'))
    throw new Error('Clock model assets are missing their index')
  if (clocks) {
    const clockAsset = manifest.assets.find((asset) => asset.path === 'clock-data.json')
    const referencedModels = new Set()
    const clockIds = new Set()
    if (clockAsset?.kind !== 'clock') throw new Error('Invalid clock index asset')
    for (const cell of clocks.cellTypes) {
      if (clockIds.has(cell.id) || (cell.scope === 'cross_cell_type'
        ? cell.id !== 'cross-cell-type-mean'
        : !de.cellTypes.some((item) => item.id === cell.id && item.label === cell.label)))
        throw new Error('Clock results reference an unknown cell type')
      clockIds.add(cell.id)
      const methodIds = new Set()
      for (const method of cell.methods) {
        const expectedPath = clockModelAssetPath(cell.id, method.id)
        const modelAsset = manifest.assets.find((asset) => asset.path === method.asset)
        const model = parsed.get(method.asset)
        if (methodIds.has(method.id) || method.asset !== expectedPath || modelAsset?.kind !== 'clockModel' ||
            !model || model.cellTypeId !== cell.id || model.methodId !== method.id)
          throw new Error(`Clock results reference an invalid model: ${cell.id}/${method.id}`)
        methodIds.add(method.id)
        referencedModels.add(method.asset)
      }
    }
    if (manifest.assets.some((asset) => asset.kind === 'clockModel' && !referencedModels.has(asset.path)))
      throw new Error('Unreferenced clock model asset')
  }
  const labels = new Set(),
    ids = new Set()
  for (const cell of de.cellTypes) {
    if (ids.has(cell.id) || labels.has(cell.label)) throw new Error('Duplicate cell type ID or label')
    ids.add(cell.id)
    labels.add(cell.label)
    if (cell.previewAsset && !seen.has(cell.previewAsset))
      throw new Error(`Missing differential-expression preview: ${cell.previewAsset}`)
    if (new Set(cell.models.map((model) => model.id)).size !== cell.models.length)
      throw new Error(`Duplicate model ID: ${cell.id}`)
    for (const model of cell.models) {
      if (model.downloadAsset && (model.model === 'per_study' || !seen.has(model.downloadAsset)))
        throw new Error(`Invalid combined-model download reference: ${cell.id}/${model.id}`)
      if (!model.recordsAsset) continue
      const chunk = parsed.get(model.recordsAsset)
      if (
        !chunk ||
        chunk.cellTypes.length !== 1 ||
        chunk.cellTypes[0].id !== cell.id ||
        chunk.cellTypes[0].models.length !== 1 ||
        chunk.cellTypes[0].models[0].id !== model.id ||
        model.records.length !== 0
      )
        throw new Error(`Invalid DE model reference: ${cell.id}/${model.id}`)
      if (
        model.recordBuckets &&
        chunk.cellTypes[0].models[0].records.some((record) => record.padj === null || record.padj >= 0.05)
      )
        throw new Error(`Unfiltered default DE model: ${cell.id}/${model.id}`)
      for (const target of Object.values(model.recordBuckets || {})) {
        const extra = parsed.get(target)
        if (
          !extra ||
          extra.cellTypes.length !== 1 ||
          extra.cellTypes[0].id !== cell.id ||
          extra.cellTypes[0].models.length !== 1 ||
          extra.cellTypes[0].models[0].id !== model.id
        )
          throw new Error(`Invalid DE gene chunk: ${cell.id}/${model.id}`)
      }
    }
    if (cell.trajectoryAsset && !parsed.has(cell.trajectoryAsset))
      throw new Error(`Missing age preview: ${cell.trajectoryAsset}`)
    if (cell.referenceExpressionAsset) {
      const expression = parsed.get(cell.referenceExpressionAsset)
      if (!expression || expression.cellType !== cell.label)
        throw new Error(`Invalid reference-expression asset reference: ${cell.id}`)
    }
  }
  const trajectory = parsed.get('pipeline-trajectory-index.json')
  for (const target of Object.values(trajectory?.geneBuckets || {})) {
    if (!parsed.has(target)) throw new Error('Missing trajectory gene chunk')
  }
  const trajectoryMap = parsed.get('pipeline-trajectory-map.json')
  for (const target of Object.values(trajectoryMap?.geneBuckets || {})) {
    const chunk = parsed.get(target)
    if (!chunk || chunk.release !== trajectoryMap.release)
      throw new Error(`Invalid trajectory detail reference: ${target}`)
  }
  for (const [path, preview] of parsed) {
    if (!path.startsWith('age-previews/') || !preview.geneBuckets) continue
    for (const target of Object.values(preview.geneBuckets)) {
      const chunk = parsed.get(target)
      if (!chunk || chunk.referenceBin !== preview.referenceBin || chunk.design !== preview.design)
        throw new Error(`Invalid age-preview gene chunk: ${path}`)
    }
  }
  if (trajectory && trajectory.records.some((record) => !labels.has(record.cellType)))
    throw new Error('Trajectory references an unknown cell type')
  const composition = parsed.get('pipeline-composition-index.json')
  if (
    composition &&
    [
      ...composition.records,
      ...composition.coefficients,
      ...composition.ageEstimates,
      ...composition.ageStatus,
    ].some((record) => !labels.has(record.cellType))
  )
    throw new Error('Composition references an unknown cell type')
  console.log(`Validated release ${manifest.releaseId}: ${seen.size} assets`)
  return manifest
}
if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href)
  await verifyRelease(process.argv[2] || 'site/public/data')
