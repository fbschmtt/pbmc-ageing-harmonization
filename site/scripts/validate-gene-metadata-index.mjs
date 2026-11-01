import { readFile } from 'node:fs/promises'
import { resolve } from 'node:path'
import { fileURLToPath } from 'node:url'
import { pantherFilename, pantherRelease, pantherSha256 } from './panther-release.mjs'

const root = resolve(fileURLToPath(new URL('..', import.meta.url)))
const metadata = JSON.parse(await readFile(resolve(root, 'public/data/gene-metadata.json'), 'utf8'))
const atlas = JSON.parse(await readFile(resolve(root, 'public/data/pipeline-trajectory-index.json'), 'utf8'))

function assert(condition, message) {
  if (!condition) throw new Error(`Gene metadata validation failed: ${message}`)
}

assert(metadata.release === `PANTHER ${pantherRelease}`, 'the PANTHER release must be pinned')
assert(metadata.source === pantherFilename, 'the source filename must be recorded')
assert(metadata.sourceSha256 === pantherSha256, 'the source checksum must be recorded')
assert(typeof metadata.sourceUrl === 'string', 'the source URL must be recorded')
assert(metadata.genes && typeof metadata.genes === 'object', 'genes must be an object')

const atlasGenes = new Set(atlas.genes)
for (const [gene, entry] of Object.entries(metadata.genes)) {
  assert(atlasGenes.has(gene), `${gene}: metadata must belong to a retained atlas gene`)
  for (const field of ['families', 'proteinClasses']) {
    assert(Array.isArray(entry[field]), `${gene}: ${field} must be an array`)
    assert(
      entry[field].every((term) => typeof term.id === 'string' && typeof term.label === 'string'),
      `${gene}: ${field} terms must have string IDs and labels`,
    )
    assert(
      new Set(entry[field].map((term) => term.id)).size === entry[field].length,
      `${gene}: ${field} IDs must be unique`,
    )
  }
  assert(entry.goTerms && typeof entry.goTerms === 'object', `${gene}: goTerms must be an object`)
  for (const [field, limit] of [
    ['biologicalProcesses', 3],
    ['molecularFunctions', 2],
    ['cellularComponents', 2],
  ]) {
    const terms = entry.goTerms[field]
    assert(Array.isArray(terms), `${gene}: goTerms.${field} must be an array`)
    assert(terms.length <= limit, `${gene}: goTerms.${field} must contain at most ${limit} terms`)
    assert(
      terms.every((term) => typeof term.id === 'string' && typeof term.label === 'string'),
      `${gene}: goTerms.${field} terms must have string IDs and labels`,
    )
  }
}

assert(Object.keys(metadata.genes).length > 0, 'at least one gene annotation is required')
console.log(`Validated PANTHER metadata for ${Object.keys(metadata.genes).length.toLocaleString()} genes`)
