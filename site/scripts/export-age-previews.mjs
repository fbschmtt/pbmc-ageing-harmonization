import { mkdir, readFile, unlink, writeFile } from 'node:fs/promises'
import { basename, dirname, isAbsolute, join, relative, resolve } from 'node:path'
import { parseCsv } from './parse-csv.mjs'

const number = (value) =>
  value == null || value.trim() === '' || !Number.isFinite(Number(value)) ? null : Number(value)

export async function writeAgePreview(dir, slug, metadata, output) {
  const fit = (metadata.models || []).find(
    (model) => model.purpose === 'age_trajectory' && model.status === 'complete',
  )
  const path = fit?.results_by_covariate?.age_bin
  if (!path) return null
  const file = resolve(dir, path)
  const rel = relative(dir, file)
  if (rel.startsWith('..') || isAbsolute(rel)) throw new Error(`Unsafe trajectory result path: ${path}`)
  console.log(`[site-export] Building age preview for ${slug} from ${path}`)
  const { headers, rows } = parseCsv(await readFile(file, 'utf8'))
  for (const name of ['gene', 'age_bin', 'log2FoldChange'])
    if (!headers.includes(name)) throw new Error(`${file}: missing ${name}`)
  const col = Object.fromEntries(headers.map((name, index) => [name, index]))
  const genes = Object.create(null)
  for (const row of rows) {
    const gene = row[col.gene],
      bin = row[col.age_bin]
    if (!gene || !/^\d+-\d+$/.test(bin || '')) throw new Error(`${file}: invalid gene or age bin`)
    const values = genes[gene] || (genes[gene] = [])
    if (values.some((value) => value.bin === bin)) throw new Error(`${file}: duplicated age bin for ${gene}`)
    values.push({ bin, effect: number(row[col.log2FoldChange]) })
  }
  for (const values of Object.values(genes))
    values.sort((a, b) => Number(a.bin.split('-')[0]) - Number(b.bin.split('-')[0]))
  const assetName = `${slug}.json`
  const target = join(dirname(output), 'age-previews')
  await mkdir(target, { recursive: true })
  await writeFile(
    join(target, assetName),
    JSON.stringify({
      schemaVersion: 1,
      referenceBin: fit.reference_bin || metadata.age_trajectory?.reference_bin || null,
      design: fit.design || null,
      genes,
    }),
  )
  console.log(`[site-export] Wrote age preview for ${slug} with ${Object.keys(genes).length.toLocaleString()} genes`)
  return `age-previews/${assetName}`
}

export async function cleanAgePreviews(output, activeAssets) {
  const target = join(dirname(output), 'age-previews')
  const manifest = join(target, 'manifest.json')
  let previous = []
  try {
    previous = JSON.parse(await readFile(manifest, 'utf8'))
  } catch (error) {
    if (error.code !== 'ENOENT') throw error
  }
  if (!Array.isArray(previous)) throw new Error('Invalid age-preview manifest')
  const keep = new Set(activeAssets.map((asset) => basename(asset)))
  for (const name of previous) {
    if (
      typeof name !== 'string' ||
      basename(name) !== name ||
      !name.endsWith('.json') ||
      name === 'manifest.json'
    )
      throw new Error('Invalid age-preview manifest entry')
    if (!keep.has(name)) {
      try {
        await unlink(join(target, name))
      } catch (error) {
        if (error.code !== 'ENOENT') throw error
      }
    }
  }
  if (keep.size || previous.length) {
    await mkdir(target, { recursive: true })
    await writeFile(manifest, JSON.stringify([...keep].sort()))
  }
}
