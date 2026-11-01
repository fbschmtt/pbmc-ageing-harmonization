import { validAsset } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { Atlas, TrajectoryRecord } from './contracts.mjs'
import type { Atlas } from './contracts.mjs'

export function isAtlas(value: unknown): value is Atlas {
  return validAsset('trajectory', value)
}

export async function loadAtlas(): Promise<Atlas> {
  const payload = await fetchReleaseAsset('pipeline-trajectory-index.json')
  if (!isAtlas(payload))
    throw new AssetError('invalid', 'The pipeline trajectory asset does not match the expected schema.')
  return payload
}

export async function loadTrajectoryChunk(path: string): Promise<Atlas> {
  const payload = await fetchReleaseAsset(path)
  if (!isAtlas(payload)) throw new AssetError('invalid', 'Invalid trajectory gene chunk')
  return payload
}
