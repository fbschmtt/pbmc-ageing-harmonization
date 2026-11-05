import { test } from 'node:test'
import assert from 'node:assert/strict'
import { mkdtemp, mkdir, writeFile, readFile, rm } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { exportAtlas } from './export.mjs'

test('export preserves model context, missing estimates and test warnings', async () => {
  const root = await mkdtemp(join(tmpdir(), 'atlas-export-'))
  try {
    const dir = join(root, 'per_type', 'monocyte', 'differential_expression')
    await mkdir(dir, { recursive: true })
    await writeFile(
      join(dir, 'run_metadata.json'),
      JSON.stringify({
        cell_type: 'Monocyte',
        analysis_mode: 'test_only',
        interpretation_warning: 'Synthetic data',
        models: [
          {
            model: 'merged',
            status: 'complete',
            study: null,
            design: '~ age',
            n_samples: 10,
            study_sample_counts: { study_a: 10 },
            results_by_covariate: { age: 'merged.csv' },
          },
          {
            model: 'combined',
            purpose: 'age_trajectory',
            status: 'complete',
            design: '~ age_bin + sex',
            reference_bin: '20-30',
            results_by_covariate: { age_bin: 'age_bin.csv' },
          },
          { model: 'per_study', status: 'skipped', results_by_covariate: { age: 'missing.csv' } },
        ],
      }),
    )
    await writeFile(
      join(dir, 'merged.csv'),
      'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.1,0.01,0.001,0.02,5,per year\nGENE2,, ,,,,per year\n',
    )
    await writeFile(
      join(dir, 'age_bin.csv'),
      'gene,age_bin,log2FoldChange\nGENE1,30-40,0.2\nGENE1,20-30,0\nGENE2,20-30,\n',
    )
    const out = join(root, 'atlas.json')
    await exportAtlas(join(root, 'per_type'), out)
    const release = JSON.parse(await readFile(out, 'utf8'))
    const cell = release.cellTypes[0],
      model = cell.models[0]
    assert.equal(cell.warning, 'Synthetic data')
    assert.equal(cell.analysisMode, 'test_only')
    assert.deepEqual(model.studies, { study_a: 10 })
    assert.equal(model.design, '~ age')
    assert.equal(model.contrast, 'per year')
    assert.equal(model.records[0].effect, 0.1)
    assert.equal(model.records[1].effect, null)
    assert.equal(model.records[1].padj, null)
    assert.equal(cell.models.length, 1)
    assert.equal(cell.trajectoryAsset, 'age-previews/monocyte.json')
    const preview = JSON.parse(await readFile(join(root, 'age-previews', 'monocyte.json'), 'utf8'))
    assert.equal(preview.referenceBin, '20-30')
    assert.deepEqual(preview.genes.GENE1, [
      { bin: '20-30', effect: 0 },
      { bin: '30-40', effect: 0.2 },
    ])
    assert.equal(preview.genes.GENE2[0].effect, null)
    const metadataPath = join(dir, 'run_metadata.json')
    const metadata = JSON.parse(await readFile(metadataPath, 'utf8'))
    metadata.models.find((item) => item.purpose === 'age_trajectory').status = 'skipped'
    await writeFile(metadataPath, JSON.stringify(metadata))
    await exportAtlas(join(root, 'per_type'), out)
    await assert.rejects(readFile(join(root, 'age-previews', 'monocyte.json')), { code: 'ENOENT' })
  } finally {
    await rm(root, { recursive: true, force: true })
  }
})

test('exports supported trajectory and parent-fraction summaries without donor data', async () => {
  const root = await mkdtemp(join(tmpdir(), 'atlas-views-'))
  try {
    const input = join(root, 'analysis', 'per_type')
    const type = join(input, 'monocyte', 'differential_expression')
    const trajectories = join(root, 'analysis', 'combined', 'trajectory')
    await mkdir(type, { recursive: true })
    await mkdir(trajectories, { recursive: true })
    await writeFile(
      join(type, 'run_metadata.json'),
      JSON.stringify({
        cell_type: 'Monocyte',
        analysis_mode: 'test_only',
        interpretation_warning: 'Synthetic data',
        models: [
          {
            model: 'merged',
            status: 'complete',
            design: '~ age',
            n_samples: 10,
            results_by_covariate: { age: 'merged.csv' },
          },
        ],
      }),
    )
    await writeFile(join(type, 'merged.csv'), 'gene,log2FoldChange,padj,contrast\nGENE1,0.1,0.01,per year\n')
    await writeFile(
      join(type, 'reference_expression.csv'),
      'cell_type,gene,reference_bin,n_studies,n_samples,reference_cpm\nMonocyte,GENE1,20-30,2,11,123.4\n',
    )
    await writeFile(
      join(type, 'fraction_model_diagnostics.tsv'),
      'cell_type\taifi_l1_parent\tresponse\tstudy\tn_samples\tmean_fraction\tresidual_sd\tbinomial_sampling_sd\tcovariates\nMonocyte\tMyeloid\tfraction_within_aifi_l1_parent\tstudy_a\t10\t0.25\t0.12\t0.03\tage;sex\n',
    )
    await writeFile(
      join(type, 'fraction_model_coefficients.tsv'),
      'cell_type\taifi_l1_parent\tresponse\tstudy\tcovariate\tterm\tcontrast\testimate\tci_low\tci_high\tn_samples\nMonocyte\tMyeloid\tfraction_within_aifi_l1_parent\tstudy_a\tage\tage\tper decade\t0.02\t-0.01\t0.05\t10\n',
    )
    await writeFile(
      join(type, 'fraction_age_bins.tsv'),
      'cell_type\taifi_l1_parent\tstudy\tdenominator\tage_bin\treference_bin\tn_samples_bin\testimate\tci_low\tci_high\treference_mean_fraction\treference_adjusted_fraction\nMonocyte\tMyeloid\tstudy_a\tl1_parent\t30-40\t20-30\t8\t0.1\t0.05\t0.15\t\t\nMonocyte\tMyeloid\tstudy_a\tl1_parent\t20-30\t20-30\t8\t0\t0\t0\t0.25\t0.3\n',
    )
    await writeFile(
      join(type, 'fraction_age_bins_metadata.json'),
      JSON.stringify([
        { study: 'study_b', denominator: 'l1_parent', status: 'skipped', reason: 'insufficient support' },
      ]),
    )
    const bins = ['20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90']
    await writeFile(join(trajectories, 'analysis_metadata.json'), JSON.stringify({ common_bins: bins }))
    await writeFile(
      join(trajectories, 'cross_cell_type_trajectory_clusters.csv'),
      [
        'gene,cell_type,trajectory_id,trajectory_cluster,umap_1,umap_2,omnibus_fdr,' +
          bins.join(',') +
          ',' +
          bins.map((bin) => `raw_${bin}`).join(','),
        'GENE1,Monocyte,GENE1 | Monocyte,2,1,-1,0.001,0,1,2,3,4,5,6,0,-1,1,2,3,4,5',
        'GENE2,Unknown,GENE2 | Unknown,2,1,-1,0.001,0,1,2,3,4,5,6,0,-1,1,2,3,4,5',
      ].join('\n'),
    )
    const out = join(root, 'data', 'de-index.json')
    const counts = await exportAtlas(input, out)
    assert.equal(counts.trajectory.input, 2)
    assert.equal(counts.trajectory.exported, 1)
    assert.deepEqual(counts.trajectory.excluded, { cell_type_not_exported: 1 })
    const atlas = JSON.parse(await readFile(join(root, 'data', 'pipeline-trajectory-index.json'), 'utf8'))
    const catalogue = JSON.parse(await readFile(out, 'utf8'))
    const reference = JSON.parse(
      await readFile(join(root, 'data', catalogue.cellTypes[0].referenceExpressionAsset), 'utf8'),
    )
    assert.equal(reference.referenceBin, '20-30')
    assert.deepEqual(reference.genes.GENE1, { cpm: 123.4, studies: 2, samples: 11 })
    const preview = await readFile(join(root, 'data', 'pipeline-trajectory-preview.svg'), 'utf8')
    assert.match(preview, /<circle /)
    assert.equal(atlas.records[0].valueAtMaxAbsChange, 5)
    assert.deepEqual(atlas.records[0].trajectory, [0, 1, 2, 3, 4, 5, 6])
    assert.equal(atlas.records[0].referenceExpression, 123.4)
    assert.equal(atlas.records[0].referenceExpressionStudies, 2)
    assert.equal(atlas.records[0].referenceExpressionSamples, 11)
    assert.equal(atlas.warning, 'Synthetic data')
    const composition = JSON.parse(
      await readFile(join(root, 'data', 'pipeline-composition-index.json'), 'utf8'),
    )
    assert.deepEqual(composition.records[0], {
      cellType: 'Monocyte',
      parent: 'Myeloid',
      study: 'study_a',
      samples: 10,
      meanFraction: 0.25,
      residualSd: 0.12,
      binomialSamplingSd: 0.03,
      covariates: 'age;sex',
      warning: 'Synthetic data',
    })
    assert.deepEqual(composition.coefficients[0], {
      cellType: 'Monocyte',
      parent: 'Myeloid',
      study: 'study_a',
      covariate: 'age',
      term: 'age',
      contrast: 'per decade',
      estimate: 0.02,
      low: -0.01,
      high: 0.05,
      samples: 10,
    })
    assert.deepEqual(composition.ageEstimates[0], {
      cellType: 'Monocyte',
      parent: 'Myeloid',
      group: 'study_a',
      denominator: 'l1_parent',
      ageBin: '30-40',
      referenceBin: '20-30',
      referenceMeanFraction: null,
      referenceAdjustedFraction: null,
      estimate: 0.1,
      low: 0.05,
      high: 0.15,
      samples: 8,
    })
    assert.equal(composition.ageEstimates[1].referenceMeanFraction, 0.25)
    assert.equal(composition.ageEstimates[1].referenceAdjustedFraction, 0.3)
    assert.equal(composition.ageStatus[0].reason, 'insufficient support')
    await writeFile(join(trajectories, 'analysis_metadata.json'), JSON.stringify({ common_bins: ['20-30'] }))
    await exportAtlas(input, out)
    await assert.rejects(readFile(join(root, 'data', 'pipeline-trajectory-index.json')), { code: 'ENOENT' })
  } finally {
    await rm(root, { recursive: true, force: true })
  }
})
