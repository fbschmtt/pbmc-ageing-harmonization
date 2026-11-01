import { validAsset, type DifferentialRelease } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { Cell, Gene, Model } from './contracts.mjs'
export async function loadDifferential(): Promise<DifferentialRelease> {
  const value = await fetchReleaseAsset('de-index.json')
  if (!validAsset('de', value)) throw new AssetError('invalid', 'Invalid DE asset')
  return value
}

export async function loadModelRecords(path: string, cellId: string, modelId: string) {
  const value = await fetchReleaseAsset(path)
  if (
    !validAsset('de', value) ||
    value.cellTypes.length !== 1 ||
    value.cellTypes[0].id !== cellId ||
    value.cellTypes[0].models.length !== 1 ||
    value.cellTypes[0].models[0].id !== modelId
  )
    throw new AssetError('invalid', 'Invalid DE model asset')
  return value.cellTypes[0].models[0].records
}
