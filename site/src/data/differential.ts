import { validAsset, type DifferentialRelease } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { Cell, Gene, Model } from './contracts.mjs'
export type ReferenceExpression = {
  schemaVersion: number
  releaseId?: string
  cellType: string
  referenceBin: string
  method: string
  genes: Record<string, { cpm: number; studies: number; samples: number }>
}
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

export async function loadReferenceExpression(path: string, cellType: string): Promise<ReferenceExpression> {
  const value = await fetchReleaseAsset(path)
  if (!validAsset('referenceExpression', value) || value.cellType !== cellType)
    throw new AssetError('invalid', 'Invalid reference-expression asset')
  return value
}

export async function preloadCombinedAgeModel(cellType?: string | null): Promise<void> {
  const release = await loadDifferential()
  const cell =
    release.cellTypes.find((item) => item.label === cellType) ??
    release.cellTypes.find((item) => item.label.toLowerCase() === 'naive cd8 t cell') ??
    release.cellTypes[0]
  const model = cell?.models.find((item) => item.model === 'merged' && item.covariate === 'age')
  if (!cell || !model?.recordsAsset) return
  await Promise.allSettled([
    loadModelRecords(model.recordsAsset, cell.id, model.id),
    ...(cell.referenceExpressionAsset
      ? [loadReferenceExpression(cell.referenceExpressionAsset, cell.label)]
      : []),
  ])
}
