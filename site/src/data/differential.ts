import { validAsset, type DifferentialRelease } from './contracts.mjs'
import { AssetError, fetchReleaseAsset, prefetchReleaseAsset } from './release'
import { geneBucket } from './geneBucket.mjs'

export type { Cell, Gene, Model } from './contracts.mjs'
export type ReferenceExpression = {
  schemaVersion: number
  releaseId?: string
  cellType: string
  referenceBin: string
  method: string
  genes: Record<string, { cpm: number; studies: number; samples: number }>
}
const combinedAgeEffectCache = new Map<string, Promise<number | null>>()

export function loadCombinedAgeEffect(cellType: string, gene: string): Promise<number | null> {
  const key = `${cellType.toLocaleLowerCase()}\0${gene.toLocaleUpperCase()}`
  const cached = combinedAgeEffectCache.get(key)
  if (cached) return cached
  const request = (async () => {
    const release = await loadDifferential()
    const cell = release.cellTypes.find((item) => item.label === cellType)
    const model = cell?.models.find((item) => item.model === 'merged' && item.covariate === 'age')
    if (!cell || !model?.recordsAsset) return null
    const findEffect = (records: Awaited<ReturnType<typeof loadModelRecords>>) =>
      records.find((record) => record.gene.toLocaleUpperCase() === gene.toLocaleUpperCase())?.effect ?? null
    const records = await loadModelRecords(model.recordsAsset, cell.id, model.id)
    const direct = findEffect(records)
    if (direct !== null) return direct
    const bucket = model.recordBuckets?.[geneBucket(gene)]
    return bucket ? findEffect(await loadModelRecords(bucket, cell.id, model.id)) : null
  })()
  combinedAgeEffectCache.set(key, request)
  request.catch(() => combinedAgeEffectCache.delete(key))
  return request
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
  const preferred =
    release.cellTypes.find((item) => item.label === cellType) ??
    release.cellTypes.find((item) => item.label.toLowerCase() === 'naive cd8 t cell') ??
    release.cellTypes[0]
  const combinedAgeModels = release.cellTypes.flatMap((cell) => {
    const model = cell.models.find((item) => item.model === 'merged' && item.covariate === 'age')
    return model?.recordsAsset ? [{ cell, model }] : []
  })
  if (!combinedAgeModels.length) return
  const ranked = [...combinedAgeModels].sort(
    (a, b) =>
      Object.keys(b.model.studies).length - Object.keys(a.model.studies).length ||
      b.model.samples - a.model.samples,
  )
  const first = combinedAgeModels.find(({ cell }) => cell === preferred) ?? ranked[0]
  const warm = [first, ...ranked.filter((item) => item !== first).slice(0, 5)]
  await Promise.allSettled(
    warm.map(({ cell, model }, index) =>
      index === 0
        ? loadModelRecords(model.recordsAsset!, cell.id, model.id)
        : prefetchReleaseAsset(model.recordsAsset!),
    ),
  )
}
