import { validAsset, type CrossTypeIndex, type CrossTypeGene, type CrossTypeGenes } from './contracts.mjs'
import { geneBucket } from './geneBucket.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { CrossTypeIndex, CrossTypeGene }
export async function loadCrossType(): Promise<CrossTypeIndex> {
  const data = await fetchReleaseAsset('cross-type-index.json')
  if (!validAsset('crossType', data)) throw new AssetError('invalid', 'Invalid cross-type summary')
  return data
}
export async function loadCrossTypeGene(index: CrossTypeIndex, gene: string): Promise<CrossTypeGene> {
  const path = index.geneBuckets[geneBucket(gene)]
  if (!path) throw new AssetError('invalid', 'Missing cross-type gene chunk')
  const data = await fetchReleaseAsset(path)
  if (!validAsset('crossTypeGenes', data)) throw new AssetError('invalid', 'Invalid cross-type gene results')
  const detail = (data as CrossTypeGenes).genes[gene]
  if (!detail) throw new AssetError('invalid', 'Missing cross-type gene results')
  return detail
}
