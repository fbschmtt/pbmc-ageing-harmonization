import { validAsset, type ClockIndex } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export async function loadClockIndex(): Promise<ClockIndex> {
  const value = await fetchReleaseAsset('clock-data.json')
  if (!validAsset('clock', value)) throw new AssetError('invalid', 'Invalid transcriptomic clock data')
  return value
}
