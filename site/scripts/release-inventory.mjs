import { readFile, readdir } from 'node:fs/promises'
import { dirname, join, resolve } from 'node:path'
import { exists } from './release-sources.mjs'

export async function releaseInventory(input, de) {
  const missingMetadata = []
  for (const entry of await readdir(input, { withFileTypes: true })) {
    if (
      entry.isDirectory() &&
      (await exists(join(input, entry.name, 'differential_expression'))) &&
      !(await exists(join(input, entry.name, 'differential_expression/run_metadata.json')))
    )
      missingMetadata.push(entry.name)
  }
  if (missingMetadata.length)
    throw new Error(
      `Differential-expression directories missing run_metadata.json: ${missingMetadata.join(', ')}`,
    )
  const inventory = []
  for (const entry of await readdir(input, { withFileTypes: true })) {
    if (!entry.isDirectory()) continue
    const path = join(input, entry.name, 'differential_expression/run_metadata.json')
    if (!(await exists(path))) continue
    const meta = JSON.parse(await readFile(path, 'utf8'))
    const cell = de.cellTypes.find((item) => item.id === entry.name)
    const artifacts = {}
    for (const name of [
      'fraction_model_diagnostics.tsv',
      'fraction_model_coefficients.tsv',
      'fraction_age_bins.tsv',
      'fraction_age_bins_metadata.json',
    ])
      artifacts[name] = (await exists(join(dirname(path), name))) ? 'present' : 'absent'
    inventory.push({
      id: entry.name,
      label: meta.cell_type,
      analysisMode: meta.analysis_mode,
      syntheticTestData: typeof meta.synthetic_test_data === 'boolean' ? meta.synthetic_test_data : null,
      warning: meta.interpretation_warning,
      exported: !!cell,
      artifacts,
      models: (meta.models || []).map((model) => ({
        model: model.model,
        study: model.study || null,
        purpose: model.purpose || null,
        status: model.status,
        reason: model.reason || null,
        covariates: Object.keys(model.results_by_covariate || {}),
      })),
    })
  }
  inventory.sort((a, b) => a.id.localeCompare(b.id))
  const deManifestPath = resolve(
    input,
    '..',
    'combined',
    'differential_expression',
    'differential_expression.json',
  )
  let analysis = null
  if (await exists(deManifestPath)) {
    const meta = JSON.parse(await readFile(deManifestPath, 'utf8'))
    const expected = Object.keys(meta.results_by_cell_type || {})
    const missing = expected.filter((id) => !inventory.some((item) => item.id === id))
    if (missing.length)
      throw new Error(`Cell types in the DE manifest lack run metadata: ${missing.join(', ')}`)
    analysis = {
      timestamp: meta.timestamp_utc || null,
      configurationSha256: meta.configuration_sha256 || null,
      method: meta.method || null,
      methodVersion: meta.method_version || null,
      splitBy: meta.split_by || null,
      analysisMode: meta.analysis_mode || null,
      syntheticTestData: typeof meta.synthetic_test_data === 'boolean' ? meta.synthetic_test_data : null,
      sampleInclusion: meta.sample_inclusion_applied || null,
      geneUniverse: meta.gene_universe || null,
      effectScales: meta.effect_scales || null,
      expectedCellTypes: expected,
    }
  }
  return { inventory, analysis }
}
