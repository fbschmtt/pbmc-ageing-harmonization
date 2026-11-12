import { afterEach, expect, it } from 'vitest'
import { mkdtemp, readFile, writeFile, rm } from 'node:fs/promises'
import { createHash } from 'node:crypto'
import { gzipSync } from 'node:zlib'
import { join } from 'node:path'
import { tmpdir } from 'node:os'
import { makeFixture } from '../tests/browser/fixture.mjs'
import { verifyRelease } from './verify-release.mjs'
import { geneBucket } from '../src/data/geneBucket.mjs'

const roots = []
afterEach(async () => {
  for (const root of roots.splice(0)) await rm(root, { recursive: true, force: true })
})
async function fixture() {
  const root = await mkdtemp(join(tmpdir(), 'cross-type-release-'))
  roots.push(root)
  const { publicDir } = await makeFixture(root)
  return join(publicDir, 'data')
}
const sha = (bytes) => createHash('sha256').update(bytes).digest('hex')
async function amend(directory, name, change) {
  const manifest = JSON.parse(await readFile(join(directory, 'release.json'), 'utf8'))
  const asset = manifest.assets.find((row) => row.path === name)
  const value = JSON.parse(await readFile(join(directory, name), 'utf8'))
  change(value)
  const bytes = Buffer.from(JSON.stringify(value))
  await writeFile(join(directory, name), bytes)
  Object.assign(asset, { bytes: bytes.length, sha256: sha(bytes) })
  if (asset.gzipPath) {
    const compressed = gzipSync(bytes)
    await writeFile(join(directory, asset.gzipPath), compressed)
    Object.assign(asset, { gzipBytes: compressed.length, gzipSha256: sha(compressed) })
  }
  const manifestBytes = Buffer.from(JSON.stringify(manifest))
  await writeFile(join(directory, 'release.json'), manifestBytes)
  await writeFile(join(directory, 'release.json.gz.bin'), gzipSync(manifestBytes))
}
it('rejects recurrence contradictions even when asset hashes are valid', async () => {
  const directory = await fixture()
  await amend(directory, 'cross-type-index.json', (index) => {
    index.genes.find((row) => row.gene === 'GENE1').significant = 1
  })
  await expect(verifyRelease(directory)).rejects.toThrow('Cross-type summary disagrees')
})
it('rejects missing gene details even when the chunk is valid and correctly hashed', async () => {
  const directory = await fixture()
  await amend(directory, `cross-type-genes/${geneBucket('GENE1')}.json`, (chunk) => {
    delete chunk.genes.GENE1
  })
  await expect(verifyRelease(directory)).rejects.toThrow('Cross-type summary disagrees')
})
