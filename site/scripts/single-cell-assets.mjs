import { createHash } from 'node:crypto'
import { readFile } from 'node:fs/promises'
import { join } from 'node:path'

export async function readSingleCellAssets(directory, cellType) {
  let manifest
  try {
    manifest = JSON.parse(await readFile(join(directory, 'site-assets.json'), 'utf8'))
  } catch (error) {
    if (error.code === 'ENOENT')
      throw new Error(`${directory}: site-assets.json is missing; rerender this cell-type report`)
    throw error
  }
  if (manifest.schemaVersion !== 1 || typeof manifest.cellType !== 'string' ||
      (cellType !== undefined && manifest.cellType !== cellType) ||
      !manifest.assets || typeof manifest.assets !== 'object' || Array.isArray(manifest.assets) ||
      !['fractionHasStudyTrends', 'fractionIncludesZeros'].every((key) => typeof manifest[key] === 'boolean') ||
      !['ageDeNote', 'fractionNote'].every((key) => manifest[key] === null || typeof manifest[key] === 'string') ||
      !(manifest.fractionDetails === null || Array.isArray(manifest.fractionDetails)))
    throw new Error(`Invalid single-cell asset manifest: ${directory}`)
  const references = [manifest.image, manifest.ageDeImage,
    ...(manifest.fractionDetails || []).map((detail) => detail.image)].filter((path) => path !== null)
  if (references.some((path) => typeof path !== 'string' || !Object.hasOwn(manifest.assets, path)))
    throw new Error(`Unindexed single-cell image: ${directory}`)
  for (const [path, record] of Object.entries(manifest.assets)) {
    if (!/^site-assets\/[a-f0-9]{64}\.png$/.test(path) ||
        typeof record?.sha256 !== 'string')
      throw new Error(`Invalid single-cell image path: ${path}`)
    const bytes = await readFile(join(directory, path))
    if (!bytes.subarray(0, 8).equals(Buffer.from('89504e470d0a1a0a', 'hex')) ||
        createHash('sha256').update(bytes).digest('hex') !== record.sha256)
      throw new Error(`Single-cell image checksum mismatch: ${path}`)
  }
  return manifest
}
