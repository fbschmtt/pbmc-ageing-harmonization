export type ReleaseManifest = {
  schemaVersion: number
  releaseId: string
  generatedAt: string
  analysis?: {
    timestamp?: string | null
    analysisMode?: string | null
    syntheticTestData?: boolean | null
    method?: string | null
    methodVersion?: string | null
    splitBy?: string | null
    configurationSha256?: string | null
    sampleInclusion?: {
      minimum_age_years_inclusive?: number
      minimum_cells_per_pseudobulk_inclusive?: number
      minimum_cells_filter_enabled?: boolean
    } | null
  } | null
  inventory?: {
    id: string
    label: string
    analysisMode?: string | null
    syntheticTestData?: boolean | null
    exported: boolean
  }[]
  provenance?: {
    status: 'consistent' | 'partial' | 'unavailable'
    checkedCellTypes: string[]
    uncheckedCellTypes?: string[]
    note: string
  }
  context?: {
    studies: {
      id: string
      publicationDoi: string | null
      annotation: string | null
      technology: string | null
      inclusion: string | null
      siteComments?: string | null
      licence: string | null
    }[]
    studyReference: { sha256: string | null; note: string }
    ageSupport: {
      commonBins: string[] | null
      referenceBinStartAge: number | null
      minimumSamplesPerBin: number | null
      cellTypes: {
        cellType: string
        status: string
        reason: string | null
        ageBinCounts: Record<string, number> | null
      }[]
    } | null
  }
  pathwaySources?: { path: string; sha256: string }[]
  assets: {
    path: string
    kind: string
    bytes: number
    sha256: string
    gzipPath?: string
    gzipBytes?: number
    gzipSha256?: string
    label?: string
    study?: string
    warning?: string | null
    unavailableCompanions?: string[]
  }[]
}
export type AssetFailureKind = 'unavailable' | 'invalid' | 'mismatch'
export class AssetError extends Error {
  constructor(
    public readonly kind: AssetFailureKind,
    message: string,
  ) {
    super(message)
    this.name = 'AssetError'
  }
}
export function assetFailureKind(error: unknown): AssetFailureKind {
  return error instanceof AssetError ? error.kind : 'unavailable'
}
async function readJson(response: Response): Promise<unknown> {
  try {
    return await response.json()
  } catch {
    throw new AssetError('invalid', 'Invalid JSON asset')
  }
}
const object = (value: unknown): value is Record<string, unknown> =>
  typeof value === 'object' && value !== null && !Array.isArray(value)
const nullableString = (value: unknown) => value === null || typeof value === 'string'
const nullableNumber = (value: unknown) =>
  value === null || (typeof value === 'number' && Number.isFinite(value))
function validContext(value: unknown): boolean {
  if (!object(value) || !Array.isArray(value.studies) || !object(value.studyReference)) return false
  if (!nullableString(value.studyReference.sha256) || typeof value.studyReference.note !== 'string')
    return false
  if (
    !value.studies.every(
      (study: unknown) =>
        object(study) &&
        typeof study.id === 'string' &&
        (study.siteComments === undefined || nullableString(study.siteComments)) &&
        ['publicationDoi', 'annotation', 'technology', 'inclusion', 'licence'].every((key) =>
          nullableString(study[key]),
        ),
    )
  )
    return false
  const age = value.ageSupport
  if (age === null) return true
  if (
    !object(age) ||
    !(
      age.commonBins === null ||
      (Array.isArray(age.commonBins) && age.commonBins.every((bin: unknown) => typeof bin === 'string'))
    ) ||
    !nullableNumber(age.referenceBinStartAge) ||
    !nullableNumber(age.minimumSamplesPerBin) ||
    !Array.isArray(age.cellTypes)
  )
    return false
  return age.cellTypes.every(
    (cell: unknown) =>
      object(cell) &&
      typeof cell.cellType === 'string' &&
      typeof cell.status === 'string' &&
      nullableString(cell.reason) &&
      (cell.ageBinCounts === null ||
        (object(cell.ageBinCounts) &&
          Object.values(cell.ageBinCounts).every(
            (count) => typeof count === 'number' && Number.isFinite(count),
          ))),
  )
}
function validProvenance(value: unknown): boolean {
  return (
    object(value) &&
    ['consistent', 'partial', 'unavailable'].includes(String(value.status)) &&
    typeof value.note === 'string' &&
    Array.isArray(value.checkedCellTypes) &&
    value.checkedCellTypes.every((id: unknown) => typeof id === 'string') &&
    (value.uncheckedCellTypes === undefined ||
      (Array.isArray(value.uncheckedCellTypes) &&
        value.uncheckedCellTypes.every((id: unknown) => typeof id === 'string')))
  )
}
let pending: Promise<ReleaseManifest> | undefined

async function fetchManifest(): Promise<unknown> {
  if (typeof DecompressionStream !== 'undefined') {
    try {
      const response = await fetch(`${import.meta.env.BASE_URL}data/release.json.gz.bin`, {
        cache: 'no-cache',
      })
      if (response.ok && response.body)
        return JSON.parse(
          await new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).text(),
        )
    } catch {
      /* Older exports do not have a compressed manifest. */
    }
  }
  const response = await fetch(`${import.meta.env.BASE_URL}data/release.json`, { cache: 'no-cache' })
  if (!response.ok) throw new AssetError('unavailable', 'Release manifest unavailable')
  return readJson(response)
}

export async function loadReleaseManifest(): Promise<ReleaseManifest> {
  if (!pending)
    pending = fetchManifest()
      .then(async (value) => {
        if (typeof value !== 'object' || value === null)
          throw new AssetError('invalid', 'Invalid release manifest')
        const data = value as ReleaseManifest
        if (
          data.schemaVersion !== 1 ||
          typeof data.releaseId !== 'string' ||
          typeof data.generatedAt !== 'string' ||
          (data.context !== undefined && !validContext(data.context)) ||
          (data.pathwaySources !== undefined &&
            (!Array.isArray(data.pathwaySources) ||
              !data.pathwaySources.every(
                (item) => item && typeof item.path === 'string' && typeof item.sha256 === 'string',
              ))) ||
          (data.provenance !== undefined && !validProvenance(data.provenance)) ||
          (data.inventory !== undefined &&
            (!Array.isArray(data.inventory) ||
              !data.inventory.every(
                (cell) =>
                  cell &&
                  typeof cell.id === 'string' &&
                  typeof cell.label === 'string' &&
                  typeof cell.exported === 'boolean' &&
                  (cell.analysisMode == null || typeof cell.analysisMode === 'string') &&
                  (cell.syntheticTestData == null || typeof cell.syntheticTestData === 'boolean'),
              ))) ||
          !Array.isArray(data.assets) ||
          !data.assets.every(
            (asset) =>
              asset &&
              typeof asset.path === 'string' &&
              typeof asset.kind === 'string' &&
              typeof asset.bytes === 'number' &&
              typeof asset.sha256 === 'string' &&
              (asset.gzipPath === undefined ||
                (asset.gzipPath === `${asset.path}.gz.bin` &&
                  typeof asset.gzipBytes === 'number' &&
                  typeof asset.gzipSha256 === 'string')) &&
              (asset.label === undefined || typeof asset.label === 'string') &&
              (asset.study === undefined || typeof asset.study === 'string') &&
              (asset.warning === undefined || asset.warning === null || typeof asset.warning === 'string') &&
              (asset.unavailableCompanions === undefined ||
                (Array.isArray(asset.unavailableCompanions) &&
                  asset.unavailableCompanions.every((item) => typeof item === 'string'))),
          )
        )
          throw new AssetError('invalid', 'Invalid release manifest')
        return data
      })
      .catch((error) => {
        pending = undefined
        throw error
      })
  return pending
}

const assetCache = new Map<string, Promise<unknown>>()
const maxCachedAssets = 48

export function fetchReleaseAsset(path: string): Promise<unknown> {
  const cached = assetCache.get(path)
  if (cached) {
    assetCache.delete(path)
    assetCache.set(path, cached)
    return cached
  }
  const request = fetchUncachedReleaseAsset(path).catch((error) => {
    if (assetCache.get(path) === request) assetCache.delete(path)
    throw error
  })
  assetCache.set(path, request)
  while (assetCache.size > maxCachedAssets) assetCache.delete(assetCache.keys().next().value!)
  return request
}

async function fetchUncachedReleaseAsset(path: string): Promise<unknown> {
  const release = await loadReleaseManifest()
  const asset = release.assets.find((item) => item.path === path)
  if (!asset) throw new AssetError('unavailable', `Asset not included in release: ${path}`)
  const compressed = asset.gzipPath && typeof DecompressionStream !== 'undefined'
  const response = await fetch(`${import.meta.env.BASE_URL}data/${compressed ? asset.gzipPath : path}`)
  if (!response.ok) throw new AssetError('unavailable', `Release asset unavailable: ${path}`)
  let value: unknown
  if (compressed && response.body) {
    try {
      value = JSON.parse(
        await new Response(response.body.pipeThrough(new DecompressionStream('gzip'))).text(),
      )
    } catch {
      throw new AssetError('invalid', `Invalid compressed asset: ${path}`)
    }
  } else value = await readJson(response)
  if (
    typeof value !== 'object' ||
    value === null ||
    !('releaseId' in value) ||
    value.releaseId !== release.releaseId
  )
    throw new AssetError('mismatch', `Release asset does not match manifest: ${path}`)
  return value
}
