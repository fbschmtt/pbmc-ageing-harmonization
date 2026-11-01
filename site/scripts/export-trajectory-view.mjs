import { writeFile } from 'node:fs/promises'
import { join, dirname } from 'node:path'
import { parseCsv } from './parse-csv.mjs'
import { num, optional, remove, exclude } from './export-view-utils.mjs'

const bins = ['20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90']

export async function exportTrajectoryView(input, output, cellTypes, counts) {
  const trajectoryOut = join(dirname(output), 'pipeline-trajectory-index.json')
  const base = join(input, '..', 'combined', 'trajectory')
  const csv = await optional(join(base, 'cross_cell_type_trajectory_clusters.csv'))
  const meta = await optional(join(base, 'analysis_metadata.json'))
  const expression = new Map()
  for (const cell of cellTypes) {
    const content = await optional(join(input, cell.id, 'differential_expression', 'reference_expression.csv'))
    if (!content) continue
    const { headers, rows } = parseCsv(content)
    const required = ['gene', 'cell_type', 'reference_bin', 'reference_cpm', 'n_studies', 'n_samples']
    if (required.some((name) => !headers.includes(name))) throw new Error('Invalid reference expression columns')
    const col = Object.fromEntries(headers.map((name, index) => [name, index]))
    for (const row of rows) {
      const value = num(row[col.reference_cpm]), studies = num(row[col.n_studies]), samples = num(row[col.n_samples])
      if (row[col.cell_type] !== cell.label || row[col.reference_bin] !== bins[0] || !row[col.gene] ||
          value === null || value < 0 || !Number.isInteger(studies) || studies < 1 ||
          !Number.isInteger(samples) || samples < studies)
        throw new Error('Invalid reference expression row')
      const key = `${cell.label}\u0000${row[col.gene]}`
      if (expression.has(key)) throw new Error('Duplicate reference expression row')
      expression.set(key, { value, studies, samples })
    }
  }
  let trajectory = []
  if (csv && meta) {
    const metadata = JSON.parse(meta)
    const { headers, rows } = parseCsv(csv)
    counts.trajectory.input = rows.length
    const required = [
      'gene',
      'cell_type',
      'trajectory_id',
      'trajectory_cluster',
      'umap_1',
      'umap_2',
      'omnibus_fdr',
      ...bins,
      ...bins.map((b) => `raw_${b}`),
    ]
    if (
      JSON.stringify(metadata.common_bins) !== JSON.stringify(bins) ||
      required.some((column) => !headers.includes(column))
    ) {
      exclude(counts.trajectory, 'unsupported_bins_or_columns', rows.length)
      console.warn(
        'Trajectory export skipped: expected seven common age bins and the associated scaled and raw fields.',
      )
    } else {
      const columns = Object.fromEntries(headers.map((name, index) => [name, index]))
      const known = new Set(cellTypes.map((cell) => cell.label))
      trajectory = rows
        .map((row) => {
          const values = bins.map((b) => num(row[columns[b]]))
          const raw = bins.map((b) => num(row[columns[`raw_${b}`]]))
          const p = num(row[columns.omnibus_fdr])
          const x = num(row[columns.umap_1]),
            y = num(row[columns.umap_2])
          if (!known.has(row[columns.cell_type])) {
            exclude(counts.trajectory, 'cell_type_not_exported')
            return null
          }
          if (
            values.some((v) => v === null) ||
            raw.some((v) => v === null) ||
            x === null ||
            y === null ||
            p === null ||
            p < 0 ||
            p > 1 ||
            !row[columns.trajectory_id] ||
            !row[columns.trajectory_cluster]
          ) {
            exclude(counts.trajectory, 'missing_or_invalid_fields')
            return null
          }
          const reference = expression.get(`${row[columns.cell_type]}\u0000${row[columns.gene]}`)
          return {
            id: row[columns.trajectory_id],
            gene: row[columns.gene],
            cellType: row[columns.cell_type],
            tissue: 'pbmc',
            x,
            y,
            adjustedPValue: p,
            referenceExpression: reference?.value ?? null,
            referenceExpressionStudies: reference?.studies ?? null,
            referenceExpressionSamples: reference?.samples ?? null,
            meanCellsPerDonor: null,
            valueAtMaxAbsChange: raw.reduce((a, b) => (Math.abs(b) > Math.abs(a) ? b : a)),
            cluster: row[columns.trajectory_cluster],
            trajectory: values,
          }
        })
        .filter(Boolean)
    }
  }
  counts.trajectory.exported = trajectory.length
  if (trajectory.length) {
    await writeFile(
      trajectoryOut,
      JSON.stringify({
        schemaVersion: 1,
        release: 'Current pipeline',
        source: 'pipeline',
        warning:
          cellTypes
            .map(
              (c) =>
                c.warning ||
                (c.analysisMode === 'test_only' ? 'Test-only results are not biological evidence.' : null),
            )
            .find(Boolean) || null,
        description:
          'Only trajectories selected for clustering are shown; p-values are adjusted omnibus age-bin tests. Shapes are scaled; maximum effects are raw log₂ fold changes relative to 20–30.',
        records: trajectory,
        tissues: ['pbmc'],
        cellTypes: [...new Set(trajectory.map((r) => r.cellType))].sort(),
        genes: [...new Set(trajectory.map((r) => r.gene))].sort(),
      }),
    )
  } else await remove(trajectoryOut)
}
