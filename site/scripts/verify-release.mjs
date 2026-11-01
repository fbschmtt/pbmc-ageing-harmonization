import { createHash } from 'node:crypto'
import { readFile } from 'node:fs/promises'
import { isAbsolute, relative, resolve } from 'node:path'
import { pathToFileURL } from 'node:url'
import { gunzipSync } from 'node:zlib'
import { validAsset } from '../src/data/contracts.mjs'

export async function verifyRelease(directory) {
  const root = resolve(directory)
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
  for (const asset of manifest.assets) {
    if (typeof asset.path !== 'string' || seen.has(asset.path))
      throw new Error('Invalid or duplicate asset path')
    seen.add(asset.path)
    const path = resolve(root, asset.path),
      rel = relative(root, path)
    if (rel.startsWith('..') || isAbsolute(rel)) throw new Error('Unsafe asset path')
    const bytes = await readFile(path)
    if (bytes.length !== asset.bytes || createHash('sha256').update(bytes).digest('hex') !== asset.sha256)
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
    const value = JSON.parse(bytes.toString('utf8'))
    parsed.set(asset.path, value)
    if (value.releaseId !== manifest.releaseId || !validAsset(asset.kind, value))
      throw new Error(`Invalid release asset: ${asset.path}`)
  }
  if (!seen.has('de-index.json')) throw new Error('Missing DE asset')
  const de = parsed.get('de-index.json')
  const labels = new Set(),
    ids = new Set()
  for (const cell of de.cellTypes) {
    if (ids.has(cell.id) || labels.has(cell.label)) throw new Error('Duplicate cell type ID or label')
    ids.add(cell.id)
    labels.add(cell.label)
    if (new Set(cell.models.map((model) => model.id)).size !== cell.models.length)
      throw new Error(`Duplicate model ID: ${cell.id}`)
    for (const model of cell.models) {
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
  }
  const trajectory = parsed.get('pipeline-trajectory-index.json')
  for (const target of Object.values(trajectory?.geneBuckets || {})) {
    if (!parsed.has(target)) throw new Error('Missing trajectory gene chunk')
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
