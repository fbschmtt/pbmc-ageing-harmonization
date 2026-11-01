import { mkdir, readFile, writeFile } from 'node:fs/promises'
import { dirname, resolve } from 'node:path'
import { fileURLToPath } from 'node:url'
import { pantherFilename, pantherRelease, pantherSha256, pantherSourceUrl } from './panther-release.mjs'
import { goTermLimits, selectInformativeTerms } from './gene-metadata-terms.mjs'

const root = resolve(dirname(fileURLToPath(import.meta.url)), '..')
const dataDir = resolve(process.env.SITE_DATA_DIR || resolve(root, 'public/data'))
const atlasSource = resolve(dataDir, 'pipeline-trajectory-index.json')
const pantherSource = resolve(root, `data-inputs/reference-data/panther/${pantherRelease}/${pantherFilename}`)
const destination = resolve(dataDir, 'gene-metadata.json')

function addTerm(terms, id, label) {
  if (!id || !label) return
  terms.set(id, { id, label })
}

function parseTerms(value) {
  const terms = new Map()
  for (const item of value.split(';')) {
    const separator = item.lastIndexOf('#')
    if (separator < 1) continue
    addTerm(terms, item.slice(separator + 1), item.slice(0, separator))
  }
  return terms
}

function sortedTerms(terms) {
  return [...terms.values()].sort((left, right) => left.label.localeCompare(right.label))
}

const atlas = JSON.parse(await readFile(atlasSource, 'utf8'))
const retainedGenes = new Set(atlas.genes)
const metadata = new Map()
const panther = await readFile(pantherSource, 'utf8')

for (const line of panther.trim().split(/\r?\n/)) {
  const columns = line.split('\t')
  const gene = columns[2]
  if (!retainedGenes.has(gene)) continue

  const entry = metadata.get(gene) ?? {
    families: new Map(),
    proteinClasses: new Map(),
    molecularFunctions: new Map(),
    biologicalProcesses: new Map(),
    cellularComponents: new Map(),
  }
  const subfamilyId = columns[3]
  const familyId = subfamilyId?.split(':')[0]
  addTerm(entry.families, familyId, columns[4])
  for (const [id, term] of parseTerms(columns[6] ?? '')) entry.molecularFunctions.set(id, term)
  for (const [id, term] of parseTerms(columns[7] ?? '')) entry.biologicalProcesses.set(id, term)
  for (const [id, term] of parseTerms(columns[8] ?? '')) entry.cellularComponents.set(id, term)
  for (const [id, term] of parseTerms(columns[9] ?? '')) entry.proteinClasses.set(id, term)
  metadata.set(gene, entry)
}

const goFields = ['molecularFunctions', 'biologicalProcesses', 'cellularComponents']
const goFrequencies = Object.fromEntries(goFields.map((field) => [field, new Map()]))
for (const entry of metadata.values()) {
  for (const field of goFields) {
    for (const id of entry[field].keys()) {
      goFrequencies[field].set(id, (goFrequencies[field].get(id) ?? 0) + 1)
    }
  }
}

const genes = Object.fromEntries(
  [...metadata.entries()]
    .sort(([left], [right]) => left.localeCompare(right))
    .map(([gene, entry]) => [
      gene,
      {
        families: sortedTerms(entry.families),
        proteinClasses: sortedTerms(entry.proteinClasses),
        goTerms: {
          molecularFunctions: selectInformativeTerms(
            entry.molecularFunctions,
            goFrequencies.molecularFunctions,
            goTermLimits.molecularFunctions,
          ),
          biologicalProcesses: selectInformativeTerms(
            entry.biologicalProcesses,
            goFrequencies.biologicalProcesses,
            goTermLimits.biologicalProcesses,
          ),
          cellularComponents: selectInformativeTerms(
            entry.cellularComponents,
            goFrequencies.cellularComponents,
            goTermLimits.cellularComponents,
          ),
        },
      },
    ]),
)
const payload = {
  schemaVersion: 1,
  release: `PANTHER ${pantherRelease}`,
  source: pantherFilename,
  sourceUrl: pantherSourceUrl,
  sourceSha256: pantherSha256,
  genes,
}

await mkdir(dirname(destination), { recursive: true })
await writeFile(destination, JSON.stringify(payload))
console.log(
  `Wrote metadata for ${Object.keys(genes).length.toLocaleString()} of ${retainedGenes.size.toLocaleString()} genes to public/data/gene-metadata.json`,
)
