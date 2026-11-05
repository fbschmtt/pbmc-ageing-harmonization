import { readFile, readdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'

async function optional(path) {
  try {
    return await readFile(path, 'utf8')
  } catch (error) {
    if (error.code === 'ENOENT') return null
    throw error
  }
}

function csvRows(text) {
  const lines = text.trim().split(/\r?\n/)
  if (!lines.length || !lines[0]) return []
  const parse = (line) => {
    const fields = []
    let field = '', quoted = false
    for (let i = 0; i < line.length; i += 1) {
      const char = line[i]
      if (char === '"' && quoted && line[i + 1] === '"') { field += '"'; i += 1 }
      else if (char === '"') quoted = !quoted
      else if (char === ',' && !quoted) { fields.push(field); field = '' }
      else field += char
    }
    fields.push(field)
    return fields
  }
  const headers = parse(lines.shift())
  return lines.filter(Boolean).map((line) => Object.fromEntries(parse(line).map((value, i) => [headers[i], value])))
}

function predictionRows(rows, { includeDonorKey = false } = {}) {
  return rows.flatMap((row) => {
    const age = Number(row.age), predictedAge = Number(row.predicted_age)
    if (!row.study || !Number.isFinite(age) || !Number.isFinite(predictedAge)) return []
    const sample = row.sample_id ?? row.sample
    const identifier = row.pseudobulk_id ?? Object.values(row)[0]
    const donorKey = sample
      ? `${row.study}::${sample}`
      : String(identifier ?? '').split('::').slice(0, -1).join('::')
    return [{
          study: row.study,
          age,
          predictedAge,
          ...(includeDonorKey && donorKey ? { donorKey } : {}),
          ...(row.n_cell_types !== undefined && row.n_cell_types !== '' && Number.isInteger(Number(row.n_cell_types))
            ? { nCellTypes: Number(row.n_cell_types) }
            : {}),
        }]
  })
}

const predictions = (rows) => predictionRows(rows)

function metrics(points) {
  const n = points.length
  if (!n) return { n: 0, mae: 0, rmse: 0, r2: null }
  const mean = points.reduce((sum, point) => sum + point.age, 0) / n
  const squared = points.reduce((sum, point) => sum + (point.age - point.predictedAge) ** 2, 0)
  const absolute = points.reduce((sum, point) => sum + Math.abs(point.age - point.predictedAge), 0)
  const total = points.reduce((sum, point) => sum + (point.age - mean) ** 2, 0)
  return { n, mae: absolute / n, rmse: Math.sqrt(squared / n), r2: total > 0 ? 1 - squared / total : null }
}

function combinePredictions(cellTypes, split) {
  if (new Set(cellTypes.map((cell) => cell.clockCohort)).size > 1) return []
  const methods = new Map()
  for (const cell of cellTypes) {
    for (const method of cell.methodFiles) {
      const rows = method[split]
      if (!rows) continue
      const predictionsByDonor = predictionRows(csvRows(rows), { includeDonorKey: true })
      for (const point of predictionsByDonor) {
        if (!point.donorKey) continue
        const grouped = methods.get(method.id) ?? new Map()
        const donors = grouped.get(point.donorKey) ?? []
        donors.push({ ...point, cellType: cell.label })
        grouped.set(point.donorKey, donors)
        methods.set(method.id, grouped)
      }
    }
  }
  return [...methods].flatMap(([id, donors]) => {
    const points = [...donors.values()].flatMap((predictions) => {
      if (predictions.length < 2) return []
      const study = predictions[0].study, age = predictions[0].age
      if (predictions.some((point) => point.study !== study || Math.abs(point.age - age) > 1e-8)) return []
      return [{
        study,
        age,
        predictedAge: predictions.reduce((sum, point) => sum + point.predictedAge, 0) / predictions.length,
        nCellTypes: predictions.length,
      }]
    })
    const name = cellTypes.flatMap((cell) => cell.methodFiles).find((method) => method.id === id)?.name ?? id
    return points.length ? [{ id, name, points, metrics: metrics(points) }] : []
  })
}

const advancedPreference = [
  'age_sex_aware_study_offset_elastic_net',
  'age_aware_study_offset_elastic_net',
  'rank_elastic_net',
  'study_zscore_elastic_net',
  'equal_study_elastic_net',
]

export async function exportClocks(input, stage, de) {
  const entries = await readdir(input, { withFileTypes: true })
  const labels = new Map(de.cellTypes.map((cell) => [cell.id, cell.label]))
  const cellTypes = []
  for (const entry of entries) {
    if (!entry.isDirectory() || !labels.has(entry.name)) continue
    const clockDir = join(input, entry.name, 'ageing_clock')
    const manifestText = await optional(join(clockDir, 'age_clock.json'))
    if (!manifestText) continue
    const manifest = JSON.parse(manifestText)
    if (manifest.kind !== 'cell_type_transcriptomic_age_clock' || !Array.isArray(manifest.methods)) continue
    const clockCohort = JSON.stringify(Object.fromEntries(
      ['input', 'configuration_sha256', 'training_studies', 'external_studies', 'cross_validation', 'normalization']
        .map((key) => [key, manifest[key]]),
    ))
    const byId = new Map(manifest.methods.map((method) => [method.method, method]))
    const chosen = []
    if (byId.has('elastic_net')) chosen.push(byId.get('elastic_net'))
    const advanced = advancedPreference.map((id) => byId.get(id)).find(Boolean)
    if (advanced && advanced.method !== 'elastic_net') chosen.push(advanced)
    const methods = []
    const methodFiles = []
    for (const method of chosen) {
      const methodDir = join(clockDir, method.method)
      const [trainingText, heldOutText] = await Promise.all([
        optional(join(methodDir, 'cv_oof_predictions.csv')),
        optional(join(methodDir, 'external_predictions.csv')),
      ])
      if (!trainingText || !heldOutText) continue
      methodFiles.push({ id: method.method, name: method.method_name, training: trainingText, heldOut: heldOutText })
      const training = predictions(csvRows(trainingText)), heldOut = predictions(csvRows(heldOutText))
      if (!training.length || !heldOut.length) continue
      methods.push({
        id: method.method,
        name: method.method_name,
        description: method.method_description,
        training: { metrics: metrics(training), predictions: training },
        heldOut: { metrics: metrics(heldOut), predictions: heldOut },
      })
    }
    if (methods.length) cellTypes.push({
      id: entry.name, label: labels.get(entry.name), methods, methodFiles, clockCohort,
    })
  }
  const crossTypeSummaryText = await optional(join(input, 'age_clock_cross_type/summary.json'))
  let crossMethods = []
  if (crossTypeSummaryText) {
    const summary = JSON.parse(crossTypeSummaryText)
    const methods = []
    for (const method of summary.methods ?? []) {
      if (!method.training?.predictions_file || !method.held_out?.predictions_file) continue
      const [trainingText, heldOutText] = await Promise.all([
        optional(join(input, 'age_clock_cross_type', method.training.predictions_file)),
        optional(join(input, 'age_clock_cross_type', method.held_out.predictions_file)),
      ])
      if (!trainingText || !heldOutText) continue
      const training = predictions(csvRows(trainingText)), heldOut = predictions(csvRows(heldOutText))
      if (!training.length || !heldOut.length) continue
      methods.push({
        id: `cross_type_mean:${method.id}`,
        name: `${method.name} · mean across cell types`,
        description: `Equal-weight mean of this clock's predictions across ${(method.training.cell_types ?? []).length} cell types for each donor. This is a cross-cell-type ensemble, not an ensemble of clock methods.`,
        training: { metrics: metrics(training), predictions: training },
        heldOut: { metrics: metrics(heldOut), predictions: heldOut },
      })
    }
    if (methods.length) cellTypes.push({
      id: 'cross-cell-type-mean',
      label: 'Across cell types · mean clock',
      scope: 'cross_cell_type',
      methods,
    })
    crossMethods = methods
  }
  // Older analysis runs may contain complete per-type predictions without the
  // cross-type aggregation summary. Derive the same donor-matched mean during
  // export so available fitted clocks appear in the site's selector.
  if (!crossMethods.length) {
    const training = combinePredictions(cellTypes, 'training')
    const heldOut = new Map(combinePredictions(cellTypes, 'heldOut').map((method) => [method.id, method]))
    const methods = training.flatMap((train) => {
      const holdout = heldOut.get(train.id)
      if (!holdout) return []
      return [{
        id: `cross_type_mean:${train.id}`,
        name: `${train.name} · mean across cell types`,
        description: 'Equal-weight mean of predictions from this clock method for each donor with predictions in at least two cell types.',
        training: { metrics: train.metrics, predictions: train.points },
        heldOut: { metrics: holdout.metrics, predictions: holdout.points },
      }]
    })
    if (methods.length) cellTypes.push({
      id: 'cross-cell-type-mean',
      label: 'Across cell types · mean clock',
      scope: 'cross_cell_type',
      methods,
    })
  }
  if (!cellTypes.length) return null
  const path = 'clock-data.json'
  // Per-type prediction CSVs are retained above only to match donors while
  // constructing a cross-cell-type mean. Do not publish those raw files or
  // their donor identifiers in the browser asset.
  const publicCellTypes = cellTypes.map(({ methodFiles: _methodFiles, clockCohort: _clockCohort, ...cell }) => cell)
  await writeFile(join(stage, path), JSON.stringify({ schemaVersion: 1, cellTypes: publicCellTypes }))
  return path
}
