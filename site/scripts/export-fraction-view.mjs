import { writeFile } from 'node:fs/promises'
import { join, dirname } from 'node:path'
import { num, optional, remove, exclude } from './export-view-utils.mjs'

export async function exportFractionView(input, output, cellTypes, counts) {
  const fractionOut = join(dirname(output), 'pipeline-composition-index.json')
  const fractions = []
  const coefficients = []
  const ageEstimates = []
  const ageStatus = []
  for (const cell of cellTypes) {
    const content = await optional(
      join(input, cell.id, 'differential_expression', 'fraction_model_diagnostics.tsv'),
    )
    if (!content) continue
    const lines = content.trim().split(/\r?\n/)
    const headers = lines.shift().split('\t')
    const required = [
      'cell_type',
      'aifi_l1_parent',
      'response',
      'study',
      'n_samples',
      'mean_fraction',
      'covariates',
    ]
    if (required.some((name) => !headers.includes(name)))
      throw new Error(`Invalid fraction diagnostics for ${cell.label}`)
    counts.fractions.input += lines.length
    for (const line of lines) {
      const row = Object.fromEntries(line.split('\t').map((value, index) => [headers[index], value]))
      const mean = num(row.mean_fraction),
        n = num(row.n_samples)
      const reason =
        row.cell_type !== cell.label
          ? 'cell_type_mismatch'
          : row.response !== 'fraction_within_aifi_l1_parent'
            ? 'unsupported_response'
            : !row.study || !row.aifi_l1_parent
              ? 'missing_identity'
              : mean === null || mean < 0 || mean > 1 || n === null
                ? 'missing_or_invalid_estimate'
                : null
      if (reason) {
        exclude(counts.fractions, reason)
        continue
      }
      fractions.push({
        cellType: cell.label,
        parent: row.aifi_l1_parent,
        study: row.study,
        samples: n,
        meanFraction: mean,
        residualSd: num(row.residual_sd),
        binomialSamplingSd: num(row.binomial_sampling_sd),
        covariates: row.covariates,
        warning:
          cell.warning ||
          (cell.analysisMode === 'test_only' ? 'Test-only results are not biological evidence.' : null),
      })
    }
  }
  for (const cell of cellTypes) {
    const content = await optional(
      join(input, cell.id, 'differential_expression', 'fraction_model_coefficients.tsv'),
    )
    if (!content) continue
    const lines = content.trim().split(/\r?\n/)
    const headers = lines.shift().split('\t')
    const required = [
      'cell_type',
      'aifi_l1_parent',
      'response',
      'study',
      'covariate',
      'term',
      'contrast',
      'estimate',
      'ci_low',
      'ci_high',
      'n_samples',
    ]
    if (required.some((name) => !headers.includes(name)))
      throw new Error(`Invalid fraction coefficients for ${cell.label}`)
    counts.coefficients.input += lines.length
    for (const line of lines) {
      const row = Object.fromEntries(line.split('\t').map((value, index) => [headers[index], value]))
      const estimate = num(row.estimate),
        low = num(row.ci_low),
        high = num(row.ci_high),
        samples = num(row.n_samples)
      const reason =
        row.cell_type !== cell.label
          ? 'cell_type_mismatch'
          : row.response !== 'fraction_within_aifi_l1_parent'
            ? 'unsupported_response'
            : !row.study || !row.aifi_l1_parent || !row.covariate || !row.contrast
              ? 'missing_identity'
              : estimate === null || low === null || high === null || samples === null
                ? 'missing_or_invalid_estimate'
                : null
      if (reason) {
        exclude(counts.coefficients, reason)
        continue
      }
      coefficients.push({
        cellType: cell.label,
        parent: row.aifi_l1_parent,
        study: row.study,
        covariate: row.covariate,
        term: row.term,
        contrast: row.contrast,
        estimate,
        low,
        high,
        samples,
      })
    }
  }
  for (const cell of cellTypes) {
    const dir = join(input, cell.id, 'differential_expression')
    const content = await optional(join(dir, 'fraction_age_bins.tsv'))
    const metadata = await optional(join(dir, 'fraction_age_bins_metadata.json'))
    if (metadata) {
      for (const status of JSON.parse(metadata)) ageStatus.push({ cellType: cell.label, ...status })
    }
    if (!content) continue
    const lines = content.trim().split(/\r?\n/)
    const headers = lines.shift().split('\t')
    const required = [
      'cell_type',
      'aifi_l1_parent',
      'study',
      'denominator',
      'age_bin',
      'reference_bin',
      'n_samples_bin',
      'estimate',
      'ci_low',
      'ci_high',
    ]
    if (required.some((name) => !headers.includes(name)))
      throw new Error(`Invalid age-bin fractions for ${cell.label}`)
    counts.ageEstimates.input += lines.length
    for (const line of lines) {
      const row = Object.fromEntries(line.split('\t').map((value, index) => [headers[index], value]))
      const estimate = num(row.estimate),
        low = num(row.ci_low),
        high = num(row.ci_high),
        samples = num(row.n_samples_bin)
      const reason =
        row.cell_type !== cell.label
          ? 'cell_type_mismatch'
          : !row.study || !row.age_bin
            ? 'missing_identity'
            : !['retained_pbmc', 'l1_parent'].includes(row.denominator)
              ? 'unsupported_denominator'
              : estimate === null || low === null || high === null || samples === null
                ? 'missing_or_invalid_estimate'
                : null
      if (reason) {
        exclude(counts.ageEstimates, reason)
        continue
      }
      ageEstimates.push({
        cellType: cell.label,
        parent: row.aifi_l1_parent,
        group: row.study,
        denominator: row.denominator,
        ageBin: row.age_bin,
        referenceBin: row.reference_bin,
        referenceMeanFraction: num(row.reference_mean_fraction),
        referenceAdjustedFraction: num(row.reference_adjusted_fraction),
        estimate,
        low,
        high,
        samples,
      })
    }
  }
  if (fractions.length || coefficients.length || ageEstimates.length)
    await writeFile(
      fractionOut,
      JSON.stringify({ schemaVersion: 1, records: fractions, coefficients, ageEstimates, ageStatus }),
    )
  else await remove(fractionOut)
  counts.fractions.exported = fractions.length
  counts.coefficients.exported = coefficients.length
  counts.ageEstimates.exported = ageEstimates.length
}
