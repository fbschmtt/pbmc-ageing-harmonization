import { validAsset, type FractionIndex } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { Fraction, FractionIndex, FractionCoefficient, FractionAgeEstimate } from './contracts.mjs'
export async function loadFractions(): Promise<FractionIndex> {
  const value = await fetchReleaseAsset('pipeline-composition-index.json')
  if (!validAsset('composition', value)) throw new AssetError('invalid', 'Invalid composition data')
  if (!value.records.length) throw new AssetError('unavailable', 'No fraction diagnostics')
  return value
}
