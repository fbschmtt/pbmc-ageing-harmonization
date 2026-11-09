import { validAsset, type ClockIndex, type ClockModelAsset } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export async function loadClockIndex(): Promise<ClockIndex> {
  const value = await fetchReleaseAsset('clock-data.json')
  if (!validAsset('clock', value)) throw new AssetError('invalid', 'Invalid transcriptomic clock data')
  return value
}

export async function loadClockModel(path: string, cellTypeId: string, methodId: string): Promise<ClockModelAsset> {
  const value = await fetchReleaseAsset(path)
  if (!validAsset('clockModel', value) || value.cellTypeId !== cellTypeId || value.methodId !== methodId)
    throw new AssetError('invalid', 'Invalid transcriptomic clock model data')
  return value
}
