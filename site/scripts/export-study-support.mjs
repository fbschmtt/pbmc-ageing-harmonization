import { readFile, writeFile } from 'node:fs/promises'
import { join } from 'node:path'

// Annotate model chunks after export, one cell type at a time, without loading
// all studies or genes into the browser.
export async function addStudySupport(stage) {
  const index = JSON.parse(await readFile(join(stage, 'de-index.json'), 'utf8'))
  for (const cell of index.cellTypes) {
    const supports = new Map()
    for (const model of cell.models.filter((item) => item.model === 'per_study' && item.study && item.contrast)) {
      const key = JSON.stringify([model.covariate, model.contrast])
      if (!supports.has(key)) supports.set(key, new Map())
      const counts = supports.get(key)
      // A study can contain multiple source rows that resolve to the same
      // displayed gene symbol. Count that gene once for the study, and treat
      // it as significant if any corresponding row passes the threshold.
      const studyGenes = new Map()
      for (const path of [model.recordsAsset, ...Object.values(model.recordBuckets || {})].filter(Boolean)) {
        const chunk = JSON.parse(await readFile(join(stage, path), 'utf8'))
        for (const record of chunk.cellTypes[0].models[0].records) {
          if (record.padj == null) continue
          const gene = record.gene
          studyGenes.set(gene, (studyGenes.get(gene) || false) || record.padj < 0.05)
        }
      }
      for (const [gene, significant] of studyGenes) {
        const value = counts.get(gene) || { significant: 0, tested: 0 }
        value.tested += 1
        if (significant) value.significant += 1
        counts.set(gene, value)
      }
    }
    for (const model of cell.models) {
      if (!model.contrast) continue
      const counts = supports.get(JSON.stringify([model.covariate, model.contrast]))
      if (!counts) continue
      for (const path of [model.recordsAsset, ...Object.values(model.recordBuckets || {})].filter(Boolean)) {
        const target = join(stage, path)
        const chunk = JSON.parse(await readFile(target, 'utf8'))
        for (const record of chunk.cellTypes[0].models[0].records)
          record.studySupport = counts.get(record.gene) || null
        await writeFile(target, JSON.stringify(chunk))
      }
    }
  }
}
