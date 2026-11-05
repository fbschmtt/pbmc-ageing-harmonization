import { readFile } from 'node:fs/promises'
import { resolve } from 'node:path'
import { exists, hash } from './release-sources.mjs'
import { parseTrajectoryMetadata } from './trajectory-metadata.mjs'

function stable(value) {
  if (Array.isArray(value)) return value.map(stable)
  if (value && typeof value === 'object')
    return Object.fromEntries(
      Object.entries(value)
        .sort(([a], [b]) => a.localeCompare(b))
        .map(([key, item]) => [key, stable(item)]),
    )
  return value
}
const same = (a, b) => JSON.stringify(stable(a)) === JSON.stringify(stable(b))

// This detects metadata contradictions, not stale tables or reports: the source
// workflow does not currently record the hashes needed to prove those links.
export async function compareAnalysisMetadata(input, inventory) {
  const dir = resolve(input, '..', 'combined', 'differential_expression')
  const path = resolve(dir, 'differential_expression.json')
  const result = {
    status: 'unavailable',
    checkedCellTypes: [],
    note: 'The combined DE manifest does not provide comparable per-type model metadata.',
  }
  if (!(await exists(path))) return result
  const combined = JSON.parse(await readFile(path, 'utf8'))
  if (!Array.isArray(combined.cell_types)) return result
  const checked = [],
    incomplete = []
  for (const cell of inventory) {
    const entries = combined.cell_types.filter((item) => item.slug === cell.id)
    if (entries.length > 1) throw new Error(`Duplicate cell type in DE manifest: ${cell.id}`)
    const entry = entries[0]
    if (!entry || !Array.isArray(entry.models)) {
      incomplete.push(cell.id)
      continue
    }
    const localDir = resolve(input, cell.id, 'differential_expression')
    const local = JSON.parse(await readFile(resolve(localDir, 'run_metadata.json'), 'utf8'))
    if (
      entry.cell_type !== local.cell_type ||
      (combined.analysis_mode && local.analysis_mode && combined.analysis_mode !== local.analysis_mode) ||
      (typeof combined.synthetic_test_data === 'boolean' &&
        typeof local.synthetic_test_data === 'boolean' &&
        combined.synthetic_test_data !== local.synthetic_test_data)
    )
      throw new Error(`DE manifest metadata disagrees for ${cell.id}`)
    if (
      local.sample_inclusion_applied &&
      combined.sample_inclusion_applied &&
      !same(local.sample_inclusion_applied, combined.sample_inclusion_applied)
    )
      throw new Error(`DE manifest sample inclusion disagrees for ${cell.id}`)
    if (!Array.isArray(local.models) || entry.models.length !== local.models.length)
      throw new Error(`DE manifest model count disagrees for ${cell.id}`)
    let complete = true
    for (let index = 0; index < entry.models.length; index++) {
      const a = local.models[index],
        b = entry.models[index]
      for (const key of ['model', 'study', 'purpose', 'status', 'design', 'n_samples']) {
        if ((a[key] ?? null) !== (b[key] ?? null))
          throw new Error(`DE manifest ${key} disagrees for ${cell.id} model ${index + 1}`)
      }
      for (const key of ['covariates', 'study_sample_counts']) {
        if (a[key] != null && b[key] != null && !same(a[key], b[key]))
          throw new Error(`DE manifest ${key} disagrees for ${cell.id} model ${index + 1}`)
      }
      if (
        !a.model ||
        !a.status ||
        (a.status === 'complete' && (!a.design || a.n_samples == null || !a.results_by_covariate))
      )
        complete = false
      const aPaths = a.results_by_covariate || {},
        bPaths = b.results_by_covariate || {}
      const names = Object.keys(aPaths).sort()
      if (JSON.stringify(names) !== JSON.stringify(Object.keys(bPaths).sort()))
        throw new Error(`DE manifest result fields disagree for ${cell.id} model ${index + 1}`)
      for (const name of names)
        if (resolve(localDir, aPaths[name]) !== resolve(dir, bPaths[name]))
          throw new Error(`DE manifest result path disagrees for ${cell.id} model ${index + 1}`)
    }
    if (complete) checked.push(cell.id)
    else incomplete.push(cell.id)
  }
  return {
    status: incomplete.length ? 'partial' : 'consistent',
    checkedCellTypes: checked,
    uncheckedCellTypes: incomplete,
    note: 'Model metadata and result paths were compared; table contents and report-to-analysis lineage are not established.',
  }
}

export async function releaseContext(input, de, studiesPath, assets = []) {
  const names = new Set()
  for (const cell of de.cellTypes)
    for (const model of cell.models) {
      if (model.study) names.add(model.study)
      for (const study of Object.keys(model.studies)) names.add(study)
    }
  const modelNames = new Set(names)
  const harmonizedNames = new Set(
    assets.filter((asset) => asset.kind === 'report' && asset.study).map((asset) => asset.study),
  )
  for (const name of harmonizedNames) names.add(name)
  const singleCellProvenancePath = resolve(
    input,
    '..',
    '..',
    'harmonization',
    'provenance',
    'single_cell_merged.json',
  )
  let singleCell = null
  if (await exists(singleCellProvenancePath))
    singleCell = JSON.parse(await readFile(singleCellProvenancePath, 'utf8'))
  const subjectsByStudy = singleCell?.subjects_by_study
  const observationsByStudy = singleCell?.observations_by_study
  let references = {},
    sourceSha256 = null
  if (await exists(studiesPath)) {
    const bytes = await readFile(studiesPath)
    sourceSha256 = hash(bytes)
    references = JSON.parse(bytes.toString('utf8')).studies || {}
  }
  const studies = [...names].sort().map((id) => {
    const entry = references[id]
    const p = entry?.provenance || {}
    return {
      id,
      inModels: modelNames.has(id),
      harmonizationAvailable: harmonizedNames.has(id),
      publicationDoi: typeof entry?.publication_doi === 'string' ? entry.publication_doi : null,
      publicationTitle: typeof entry?.publication_title === 'string' ? entry.publication_title : null,
      publicationCitation: typeof entry?.publication_citation === 'string' ? entry.publication_citation : null,
      donors: Number.isSafeInteger(subjectsByStudy?.[id]) ? subjectsByStudy[id] : null,
      cells: Number.isSafeInteger(observationsByStudy?.[id]) ? observationsByStudy[id] : null,
      annotation: typeof p.celltypes === 'string' ? p.celltypes : null,
      technology: typeof p.technology === 'string' ? p.technology : null,
      inclusion: typeof p.model_inclusion === 'string' ? p.model_inclusion : null,
      siteComments: typeof entry?.site_comments === 'string' ? entry.site_comments : null,
      licence: null,
    }
  })
  let singleCellSummary = null
  if (singleCell) {
    if (
      Number.isSafeInteger(singleCell.n_observations) &&
      singleCell.n_observations >= 0 &&
      subjectsByStudy &&
      typeof subjectsByStudy === 'object' &&
      !Array.isArray(subjectsByStudy) &&
      Object.values(subjectsByStudy).every((count) => Number.isSafeInteger(count) && count >= 0)
    ) {
      singleCellSummary = {
        studies: Object.keys(subjectsByStudy).length,
        donors: Object.values(subjectsByStudy).reduce((total, count) => total + count, 0),
        cells: singleCell.n_observations,
      }
    }
  }
  const trajectory = resolve(input, '..', 'combined', 'trajectory', 'analysis_metadata.json')
  let ageSupport = null
  if (await exists(trajectory)) {
    const meta = parseTrajectoryMetadata(await readFile(trajectory, 'utf8'))
    ageSupport = {
      commonBins: meta.common_bins ?? null,
      referenceBinStartAge: meta.reference_bin_start_age ?? null,
      minimumSamplesPerBin: meta.minimum_samples_per_bin ?? null,
      cellTypes: Array.isArray(meta.cell_type_support)
        ? meta.cell_type_support.map((item) => ({
            cellType: item.cell_type,
            status: item.status,
            reason: item.reason ?? null,
            ageBinCounts: item.age_bin_counts ?? null,
          }))
        : [],
    }
  }
  return {
    studies,
    singleCellSummary,
    studyReference: {
      sha256: sourceSha256,
      note: 'Reference metadata comes from the exporter checkout study configuration, not a verified analysis-time snapshot. Licences require review.',
    },
    ageSupport,
  }
}
