import { createHash } from 'node:crypto'
import { access, readFile, readdir } from 'node:fs/promises'
import { join, relative, resolve } from 'node:path'

export const exists = async (path) =>
  access(path).then(
    () => true,
    (e) => {
      if (e.code === 'ENOENT') return false
      throw e
    },
  )
export const hash = (contents) => createHash('sha256').update(contents).digest('hex')
export async function files(dir) {
  const entries = await readdir(dir, { withFileTypes: true })
  const result = []
  for (const entry of entries.sort((a, b) => a.name.localeCompare(b.name))) {
    const path = join(dir, entry.name)
    if (entry.isDirectory()) result.push(...(await files(path)))
    else if (entry.isFile()) result.push(path)
  }
  return result
}
export async function sourceInventory(input) {
  const result = []
  for (const entry of await readdir(input, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue
    const dir = join(input, entry.name, 'differential_expression')
    if (!(await exists(dir))) continue
    const metadataPath = join(dir, 'run_metadata.json')
    if (!(await exists(metadataPath))) continue
    const metadata = JSON.parse(await readFile(metadataPath, 'utf8'))
    const names = new Set([
      'run_metadata.json',
      'fraction_model_diagnostics.tsv',
      'fraction_model_coefficients.tsv',
      'fraction_age_bins.tsv',
      'fraction_age_bins_metadata.json',
      'reference_expression.csv',
      'reference_expression_per_study.csv',
      'reference_expression_metadata.json',
    ])
    for (const model of metadata.models || [])
      if (model.status === 'complete')
        for (const path of Object.values(model.results_by_covariate || {})) names.add(path)
    for (const name of names) {
      const path = resolve(dir, name),
        rel = relative(dir, path)
      if (rel.startsWith('..') || rel.startsWith('/')) throw new Error(`Unsafe result path: ${name}`)
      if (await exists(path)) result.push({ path: relative(input, path), sha256: hash(await readFile(path)) })
    }
  }
  const combined = resolve(input, '..', 'combined', 'trajectory')
  for (const name of ['cross_cell_type_trajectory_clusters.csv', 'analysis_metadata.json']) {
    const path = join(combined, name)
    if (await exists(path))
      result.push({ path: `../combined/trajectory/${name}`, sha256: hash(await readFile(path)) })
  }
  const deManifest = resolve(
    input,
    '..',
    'combined',
    'differential_expression',
    'differential_expression.json',
  )
  if (await exists(deManifest))
    result.push({
      path: '../combined/differential_expression/differential_expression.json',
      sha256: hash(await readFile(deManifest)),
    })
  return result.sort((a, b) => a.path.localeCompare(b.path))
}
