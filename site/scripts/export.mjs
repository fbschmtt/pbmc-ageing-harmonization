import { readdir, readFile, mkdir, writeFile } from 'node:fs/promises'
import { resolve, join, relative, isAbsolute } from 'node:path'
import { pathToFileURL } from 'node:url'
import { parseCsv } from './parse-csv.mjs'
import { exportViews } from './export-views.mjs'
import { cleanAgePreviews, writeAgePreview } from './export-age-previews.mjs'

export async function exportAtlas(inputPath, outputPath, options = {}) {
  const input = resolve(inputPath)
  const output = resolve(outputPath)
  const numeric = (value) =>
    value == null || value.trim() === '' || !Number.isFinite(Number(value)) ? null : Number(value)
  const cellTypes = []
  for (const entry of (await readdir(input, { withFileTypes: true }))
    .filter((e) => e.isDirectory())
    .sort((a, b) => a.name.localeCompare(b.name))) {
    const dir = join(input, entry.name, 'differential_expression')
    let metadata
    try {
      metadata = JSON.parse(await readFile(join(dir, 'run_metadata.json'), 'utf8'))
    } catch (error) {
      if (error.code === 'ENOENT') continue
      throw error
    }
    const models = []
    const cell = {
      id: entry.name,
      label: metadata.cell_type,
      tissue: 'pbmc',
      analysisMode: metadata.analysis_mode,
      warning: metadata.interpretation_warning,
      trajectoryAsset: null,
      models,
    }
    for (const model of metadata.models || []) {
      if (model.status !== 'complete' || model.purpose === 'age_trajectory') continue
      for (const [covariate, path] of Object.entries(model.results_by_covariate || {})) {
        const file = resolve(dir, path)
        const rel = relative(dir, file)
        if (rel.startsWith('..') || isAbsolute(rel)) throw new Error(`Unsafe result path: ${path}`)
        const { headers, rows } = parseCsv(await readFile(file, 'utf8'))
        for (const name of ['gene', 'log2FoldChange', 'padj'])
          if (!headers.includes(name)) throw new Error(`${file}: missing ${name}`)
        const col = Object.fromEntries(headers.map((name, index) => [name, index]))
        const records = rows.map((row) => ({
          gene: row[col.gene],
          effect: numeric(row[col.log2FoldChange]),
          padj: numeric(row[col.padj]),
          pvalue: numeric(row[col.pvalue]),
          se: numeric(row[col.lfcSE]),
          baseMean: numeric(row[col.baseMean]),
        }))
        if (records.some((r) => !r.gene || (r.padj !== null && (r.padj < 0 || r.padj > 1))))
          throw new Error(`${file}: invalid gene or adjusted p-value`)
        const first = rows[0] || []
        const exportedModel = {
          id: `${model.model}:${model.study || 'all'}:${covariate}`,
          model: model.model,
          study: model.study || null,
          covariate,
          contrast: first[col.contrast] || null,
          design: model.design,
          samples: model.n_samples,
          studies: model.study_sample_counts || {},
          records,
        }
        if (options.writeModel) {
          if (!models.length) cell.trajectoryAsset = await writeAgePreview(dir, entry.name, metadata, output)
          await options.writeModel(cell, exportedModel, cellTypes.length, models.length)
        }
        models.push(exportedModel)
      }
    }
    if (models.length) {
      if (!options.writeModel) cell.trajectoryAsset = await writeAgePreview(dir, entry.name, metadata, output)
      cellTypes.push(cell)
    }
  }
  if (!cellTypes.length) throw new Error(`No completed differential-expression results found in ${input}`)
  await mkdir(resolve(output, '..'), { recursive: true })
  await writeFile(
    output,
    JSON.stringify({ schemaVersion: 1, source: 'pseudobulk differential expression', cellTypes }),
  )
  await cleanAgePreviews(output, cellTypes.map((cell) => cell.trajectoryAsset).filter(Boolean))
  const counts = await exportViews(input, output, cellTypes)
  console.log(`Exported ${cellTypes.length} cell types to ${output}`)
  return counts
}
if (process.argv[1] && import.meta.url === pathToFileURL(resolve(process.argv[1])).href) {
  await exportAtlas(process.argv[2] || 'output/analysis/per_type', process.argv[3] || 'site/data/atlas.json')
}
