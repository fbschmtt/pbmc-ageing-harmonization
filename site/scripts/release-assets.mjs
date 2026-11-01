import { readFile, writeFile } from 'node:fs/promises'
import { relative } from 'node:path'
import { validAsset } from '../src/data/contracts.mjs'
import { files, hash } from './release-sources.mjs'

const kindOf = (name) =>
  name === 'de-index.json'
    ? 'de'
    : name.startsWith('de-models/') && name.endsWith('.json')
      ? 'de'
      : name.startsWith('trajectory-extras/') && name.endsWith('.json')
        ? 'trajectory'
        : name === 'pipeline-trajectory-index.json'
          ? 'trajectory'
          : name === 'pipeline-composition-index.json'
            ? 'composition'
            : name === 'gene-metadata.json'
              ? 'metadata'
              : name === 'pathway-index.json'
                ? 'pathways'
                : name.startsWith('age-previews/') && name !== 'age-previews/manifest.json'
                  ? 'preview'
                  : null

export async function stampAssets(stage, id) {
  const assets = []
  for (const path of await files(stage)) {
    const name = relative(stage, path).replaceAll('\\', '/')
    const kind = kindOf(name)
    if (!kind) continue // internal age-preview cleanup manifest
    const value = JSON.parse(await readFile(path, 'utf8'))
    if (!validAsset(kind, value)) throw new Error(`Invalid ${kind} asset: ${name}`)
    value.releaseId = id
    const contents = JSON.stringify(value)
    await writeFile(path, contents)
    assets.push({ path: name, kind, bytes: Buffer.byteLength(contents), sha256: hash(contents) })
  }
  return assets
}
