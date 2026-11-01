import { mkdir, readdir, readFile, writeFile } from 'node:fs/promises'
import { dirname, extname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const sourceDirectory = resolve(process.env.SITE_PATHWAY_DIR || resolve(root, 'public/pathways'))
const destination = resolve(process.env.SITE_DATA_DIR || resolve(root, 'public/data'), 'pathway-index.json')

let files = []
try {
  files = (await readdir(sourceDirectory)).filter((file) => extname(file).toLowerCase() === '.gmt').sort()
} catch (error) {
  if (error?.code !== 'ENOENT') throw error
}

let labels = {}
try {
  const labelsPath = resolve(sourceDirectory, 'pathway-labels.json')
  const parsed = JSON.parse(await readFile(labelsPath, 'utf8'))
  if (
    typeof parsed !== 'object' ||
    parsed === null ||
    Array.isArray(parsed) ||
    !Object.values(parsed).every((label) => typeof label === 'string')
  )
    throw new Error('pathway-labels.json must be an object of pathway IDs to labels')
  labels = parsed
} catch (error) {
  if (error?.code !== 'ENOENT') throw error
}

function displayName(id) {
  return labels[id] ?? id.replaceAll('_', ' ').replaceAll('-', ' ')
}

const pathwaySets = []
const ids = new Set()
for (const file of files) {
  const contents = await readFile(resolve(sourceDirectory, file), 'utf8')
  for (const [lineNumber, line] of contents.split(/\r?\n/).entries()) {
    if (!line.trim() || line.startsWith('#')) continue
    const [id, description, ...genes] = line.split('\t')
    if (!id || !description || genes.length === 0)
      throw new Error(`${file}:${lineNumber + 1} must contain an ID, description, and at least one gene`)
    if (ids.has(id)) throw new Error(`Duplicate gene-set ID: ${id}`)
    ids.add(id)
    pathwaySets.push({
      id,
      label: displayName(id),
      description,
      genes: [...new Set(genes.filter(Boolean).map((gene) => gene.toUpperCase()))],
      sourceFile: file,
    })
  }
}

await mkdir(dirname(destination), { recursive: true })
await writeFile(destination, JSON.stringify({ schemaVersion: 1, pathwaySets }))
console.log(`Wrote ${pathwaySets.length} pathway gene sets from ${files.length} GMT files`)
