import { readFile, writeFile } from 'node:fs/promises'
import { gzipSync } from 'node:zlib'
import { parseCsv } from './parse-csv.mjs'

const downloadColumns = [
  'gene',
  'baseMean',
  'log2FoldChange',
  'lfcSE',
  'stat',
  'pvalue',
  'padj',
]

const metadataColumns = [
  ['cell_type', 'Cell type'],
  ['model', 'Model'],
  ['study', 'Study scope'],
  ['analysis_mode', 'Analysis mode'],
  ['covariate', 'Covariate'],
  ['contrast', 'Contrast'],
  ['design', 'Design'],
  ['n_samples_model', 'Samples in model'],
  ['n_samples_before_complete_case_filter', 'Samples before complete-case filtering'],
  ['n_samples_excluded_missing_design_covariates', 'Samples excluded for missing design covariates'],
  ['study_sample_counts', 'Study sample counts'],
  ['study_site_sample_counts', 'Study-site sample counts'],
  ['interpretation_warning', 'Interpretation'],
]

const encodeField = (value) => {
  const text = value ?? ''
  return /[",\r\n]/.test(text) ? `"${text.replaceAll('"', '""')}"` : text
}

export async function writeDifferentialDownload(source, target) {
  const { headers, rows } = parseCsv(await readFile(source, 'utf8'))
  for (const required of ['gene', 'log2FoldChange', 'padj']) {
    if (!headers.includes(required)) throw new Error(`${source}: missing ${required} download column`)
  }
  const selected = downloadColumns.filter((column) => headers.includes(column))
  const indexes = selected.map((column) => headers.indexOf(column))
  const metadata = metadataColumns.flatMap(([column, label]) => {
    const index = headers.indexOf(column)
    if (index < 0 || !rows.length) return []
    const values = [...new Set(rows.map((row) => row[index]).filter(Boolean))]
    if (!values.length) return []
    if (values.length > 1) throw new Error(`${source}: ${column} varies within a combined download`)
    let value = values[0]
    if (column === 'study_sample_counts' || column === 'study_site_sample_counts') {
      try {
        const counts = JSON.parse(value)
        if (counts && typeof counts === 'object' && !Array.isArray(counts))
          value = Object.entries(counts).map(([name, count]) => `${name}: ${count}`).join('; ')
      } catch {
        // Preserve the source value if legacy output is not valid JSON.
      }
    }
    return [`# ${label}: ${value.replace(/[\r\n]+/g, ' ')}`]
  })
  const output = [
    ...metadata,
    selected,
    ...rows.map((row) => indexes.map((index) => row[index])),
  ].map((row) => Array.isArray(row) ? row.map(encodeField).join(',') : row).join('\n') + '\n'
  await writeFile(target, target.endsWith('.gz') ? gzipSync(output, { level: 9 }) : output)
}
