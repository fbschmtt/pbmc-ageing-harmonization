import { createHash } from 'node:crypto'
import { mkdir, readFile, rename, writeFile } from 'node:fs/promises'
import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'
import { pantherFilename, pantherRelease, pantherSha256, pantherSourceUrl } from './panther-release.mjs'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const destination = resolve(root, `data-inputs/reference-data/panther/${pantherRelease}/${pantherFilename}`)
const temporaryDestination = `${destination}.download`

function checksum(data) {
  return createHash('sha256').update(data).digest('hex')
}

try {
  const existing = await readFile(destination)
  if (checksum(existing) !== pantherSha256)
    throw new Error(`${destination} exists but does not match the pinned SHA-256 checksum`)
  console.log(`PANTHER ${pantherRelease} human classification is already cached and verified`)
  process.exit(0)
} catch (error) {
  if (error.code !== 'ENOENT') throw error
}

const response = await fetch(pantherSourceUrl)
if (!response.ok) throw new Error(`PANTHER download failed: HTTP ${response.status}`)
const data = Buffer.from(await response.arrayBuffer())
if (checksum(data) !== pantherSha256) throw new Error('PANTHER download does not match the pinned checksum')

await mkdir(dirname(destination), { recursive: true })
await writeFile(temporaryDestination, data)
await rename(temporaryDestination, destination)
console.log(`Cached and verified PANTHER ${pantherRelease} human classification at ${destination}`)
