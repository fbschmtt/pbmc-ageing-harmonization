import { validAsset } from './contracts.mjs'
import { AssetError, fetchReleaseAsset } from './release'

export type { GeneAnnotationTerm, GeneMetadata, GeneMetadataIndex } from './contracts.mjs'
import type { GeneMetadataIndex } from './contracts.mjs'

export async function loadGeneMetadata(): Promise<GeneMetadataIndex> {
  const payload = await fetchReleaseAsset('gene-metadata.json')
  if (!validAsset('metadata', payload))
    throw new AssetError('invalid', 'The gene metadata asset does not match the expected schema.')
  return payload
}
