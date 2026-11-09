import { mkdtemp, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { afterEach, expect, it } from 'vitest'
import { writeDifferentialDownload } from './export-de-download.mjs'

const roots = []
afterEach(async () => { await Promise.all(roots.splice(0).map((root) => rm(root, { recursive: true, force: true }))) })

it('puts fixed model metadata in a comment preamble and drops redundant columns', async () => {
  const root = await mkdtemp(join(tmpdir(), 'de-download-')); roots.push(root)
  const source = join(root, 'merged.csv'), target = join(root, 'download.csv')
  await writeFile(source, [
    'cell_type,model,study,analysis_mode,interpretation_warning,gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj,covariate,contrast,design,n_samples_model,studies_included,study_sample_counts,study_sites_included,study_site_sample_counts',
    'Monocyte,merged,all,standard,,GENE1,12,0.4,0.1,4,0.0001,0.01,age,per year,"~ study_site + age + sex",80,"study_a;study_b","{""study_a"": 40, ""study_b"": 40}","site_a;site_b","{""site_a"": 40, ""site_b"": 40}"',
  ].join('\n'))

  await writeDifferentialDownload(source, target)

  expect(await readFile(target, 'utf8')).toBe([
    '# Cell type: Monocyte',
    '# Model: merged',
    '# Study scope: all',
    '# Analysis mode: standard',
    '# Covariate: age',
    '# Contrast: per year',
    '# Design: ~ study_site + age + sex',
    '# Samples in model: 80',
    '# Study sample counts: study_a: 40; study_b: 40',
    '# Study-site sample counts: site_a: 40; site_b: 40',
    'gene,baseMean,log2FoldChange,lfcSE,stat,pvalue,padj',
    'GENE1,12,0.4,0.1,4,0.0001,0.01',
    '',
  ].join('\n'))
})
