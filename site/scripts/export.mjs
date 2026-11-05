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
  const entries = (await readdir(input, { withFileTypes: true }))
    .filter((e) => e.isDirectory())
    .sort((a, b) => a.name.localeCompare(b.name))
  let visited = 0
  for (const entry of entries) {
    visited += 1
    console.log(`[site-export] Scanning cell type ${visited}/${entries.length}: ${entry.name}`)
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
    const referencePath = join(dir, 'reference_expression.csv')
    try {
      const { headers, rows } = parseCsv(await readFile(referencePath, 'utf8'))
      const required = ['gene', 'reference_bin', 'n_studies', 'n_samples', 'reference_cpm']
      if (required.some((name) => !headers.includes(name)))
        throw new Error(`${referencePath}: missing reference-expression columns`)
      const col = Object.fromEntries(headers.map((name, index) => [name, index]))
      const referenceBin = rows[0]?.[col.reference_bin]
      const genes = {}
      for (const row of rows) {
        const gene = row[col.gene]
        const cpm = numeric(row[col.reference_cpm])
        const studies = numeric(row[col.n_studies])
        const samples = numeric(row[col.n_samples])
        if (!gene || !referenceBin || row[col.reference_bin] !== referenceBin || cpm === null || cpm < 0 ||
            !Number.isInteger(studies) || studies < 1 || !Number.isInteger(samples) || samples < studies ||
            Object.hasOwn(genes, gene))
          throw new Error(`${referencePath}: invalid or duplicate reference-expression row`)
        genes[gene] = { cpm, studies, samples }
      }
      const methodPath = join(dir, 'reference_expression_metadata.json')
      let method = 'Equal-study mean of per-study mean log1p(CPM), back-transformed with expm1'
      try {
        const metadata = JSON.parse(await readFile(methodPath, 'utf8'))
        if (typeof metadata.method === 'string' && metadata.method) method = metadata.method
      } catch (error) {
        if (error.code !== 'ENOENT') throw error
      }
      const referenceExpressionAsset = `reference-expression/${entry.name}.json`
      const referenceExpressionDirectory = resolve(output, '..', 'reference-expression')
      await mkdir(referenceExpressionDirectory, { recursive: true })
      await writeFile(
        resolve(output, '..', referenceExpressionAsset),
        JSON.stringify({ schemaVersion: 1, cellType: metadata.cell_type, referenceBin, method, genes }),
      )
      cell.referenceExpressionAsset = referenceExpressionAsset
    } catch (error) {
      if (error.code !== 'ENOENT') throw error
    }
    for (const model of metadata.models || []) {
      if (model.status !== 'complete' || model.purpose === 'age_trajectory') continue
      for (const [covariate, path] of Object.entries(model.results_by_covariate || {})) {
        const file = resolve(dir, path)
        const rel = relative(dir, file)
        if (rel.startsWith('..') || isAbsolute(rel)) throw new Error(`Unsafe result path: ${path}`)
        console.log(`[site-export] Parsing ${entry.name}: ${model.model}/${covariate} (${path})`)
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
        console.log(
          `[site-export] Parsed ${rows.length.toLocaleString()} rows for ${entry.name}: ${model.model}/${covariate}`,
        )
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
          await options.writeModel(cell, exportedModel, cellTypes.length, models.length, file)
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
