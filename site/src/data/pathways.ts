import { validAsset } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { PathwaySet } from './contracts.mjs'
import type { PathwaySet } from './contracts.mjs'

export async function loadPathwaySets(): Promise<PathwaySet[]> {
  const payload = await fetchReleaseAsset('pathway-index.json')
  if (!validAsset('pathways', payload)) throw new AssetError('invalid', 'Invalid pathway data')
  return payload.pathwaySets
}
