import { validAsset, type Preview } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'
import { geneBucket } from './geneBucket.mjs'

export type { Bin, Preview } from './contracts.mjs'
export async function loadAgePreview(asset: string, gene?: string): Promise<Preview> {
  const value = await fetchReleaseAsset(asset)
  if (!validAsset('preview', value)) throw new AssetError('invalid', 'Invalid age-bin data')
  if (value.geneBuckets && gene) {
    const path = value.geneBuckets[geneBucket(gene)]
    if (!path) return value
    const chunk = await fetchReleaseAsset(path)
    if (
      !validAsset('preview', chunk) ||
      chunk.design !== value.design ||
      chunk.referenceBin !== value.referenceBin
    )
      throw new AssetError('invalid', 'Invalid age-bin gene chunk')
    return chunk
  }
  return value
}
