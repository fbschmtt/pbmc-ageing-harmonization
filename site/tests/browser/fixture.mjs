import { mkdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'
import { createRelease } from '../../scripts/release.mjs'

const bins = ['20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90']
const write = async (path, content) => {
  await mkdir(join(path, '..'), { recursive: true })
  await writeFile(path, content)
}

// Deliberately different values make a wrong model, gene or denominator observable.
export async function makeFixture(root) {
  const input = join(root, 'analysis/per_type')
  const writeClock = async (slug, label, externalStudy) => {
    const clockDir = join(input, slug, 'ageing_clock')
    const methods = [
      { method: 'elastic_net', method_name: 'Elastic net', method_description: 'Baseline expression model.' },
      { method: 'equal_study_elastic_net', method_name: 'Equal-study elastic net', method_description: 'Equal study weighting.' },
      { method: 'rank_elastic_net', method_name: 'Rank elastic net', method_description: 'Rank-transformed expression.' },
      { method: 'study_zscore_elastic_net', method_name: 'Study z-score elastic net', method_description: 'Study-normalized expression.' },
      { method: 'age_aware_study_offset_elastic_net', method_name: 'Age-aware study-offset elastic net', method_description: 'Age-aware study adjustment.' },
    ]
    await write(join(clockDir, 'age_clock.json'), JSON.stringify({ kind: 'cell_type_transcriptomic_age_clock', cell_type: label, methods }))
    for (const method of methods) {
      const methodDir = join(clockDir, method.method)
      const offset = methods.indexOf(method) * 0.5
      const training = `pseudobulk_id,study,age,predicted_age\nsample-a,Study A,40,${42 + offset}\nsample-b,Study A,60,${58 - offset}\n`
      const external = `pseudobulk_id,study,age,predicted_age\nsample-c,${externalStudy},35,${39 + offset}\nsample-d,${externalStudy},70,${66 - offset}\n`
      await write(join(methodDir, 'cv_oof_predictions.csv'), training)
      await write(join(methodDir, 'external_predictions.csv'), external)
    }
  }
  const analysis = join(input, 'cd14-monocyte/analysis')
  await write(join(analysis, 'report.json'), JSON.stringify({ status: 'complete', cell_type: 'CD14 monocyte', n_cells: 42, n_clusters: 2, embedding: { integration: { method: 'harmony' } } }))
  await write(join(analysis, 'markers.tsv'), 'names\tscores\tlogfoldchanges\tpvals_adj\tgroup\nTESTMARKER\t3.2\t2.4\t0.01\t0\nOTHER\t-2\t-1.5\t0.1\t0\nSECOND\t3\t1\t0.01\t1\n')
  await write(join(analysis, 'executed.ipynb'), JSON.stringify({ cells: [
    { cell_type: 'markdown', source: ['### Harmony-Integrated UMAP'] },
    { cell_type: 'code', outputs: [{ data: { 'image/png': 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII=' } }] },
    { cell_type: 'code', source: ['# cluster-fraction zeros included\ncluster_fractions = sample_cluster_fractions('], outputs: [
      { data: { 'text/markdown': ['_3 sample(s) excluded from cluster fractions because the split cell-type denominator has fewer than 10 cells._'] } },
      { data: { 'text/html': ['<details><summary>Show cluster-composition curves</summary><img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII="></details>'] } },
      { data: { 'application/vnd.pbmc.cluster-fractions+json': { clusters: [
        { cluster: '0', image: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII=', fits: [{ study: null, samples: 3, meanFraction: 0.25, slopePpPerDecade: 1.5 }, { study: 'study_a', samples: 1, meanFraction: 0.5, slopePpPerDecade: null }] },
        { cluster: '1', image: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII=', fits: [{ study: null, samples: 3, meanFraction: 0.75, slopePpPerDecade: -1.5 }] },
      ] } } },
    ] },
    { cell_type: 'code', source: ['# study-specific linear trends'], outputs: [] },
  ] }))
  const dir = join(input, 'cd14-monocyte/differential_expression')
  await writeClock('cd14-monocyte', 'CD14 monocyte', 'Study B')
  const crossClock = join(input, 'age_clock_cross_type')
  await write(join(crossClock, 'summary.json'), JSON.stringify({ methods: [{
    id: 'elastic_net', name: 'Elastic net',
    training: { predictions_file: 'elastic_net/training_predictions.csv', cell_types: ['CD14 monocyte', 'Naive CD4 T cell'] },
    held_out: { predictions_file: 'elastic_net/held_out_predictions.csv' },
  }] }))
  await write(join(crossClock, 'elastic_net/training_predictions.csv'), 'donor_key,study,age,predicted_age,n_cell_types\na,Study A,40,41,2\nb,Study A,60,59,2\n')
  await write(join(crossClock, 'elastic_net/held_out_predictions.csv'), 'donor_key,study,age,predicted_age,n_cell_types\nc,Study B,35,37,2\nd,Study B,70,68,2\n')
  const trajectory = join(root, 'analysis/combined/trajectory')
  const pathwayDirectory = join(root, 'pathways')
  await mkdir(pathwayDirectory, { recursive: true })
  await write(
    join(dir, 'run_metadata.json'),
    JSON.stringify({
      cell_type: 'CD14 monocyte',
      analysis_mode: 'test_only',
      synthetic_test_data: true,
      interpretation_warning: 'Synthetic fixture: not biological evidence.',
      models: [
        {
          model: 'merged',
          status: 'complete',
          design: '~ age + sex',
          n_samples: 12,
          // The combined metadata includes 13 studies, while the fixture only
          // supplies per-gene support from three studies. The UI scale must use
          // the tested count carried by each gene, not this metadata key count.
          study_sample_counts: Object.fromEntries([
            ['study_a', 12],
            ['study_b', 12],
            ['study_c', 12],
            ...Array.from({ length: 10 }, (_, index) => [`synthetic_${index + 1}`, 12]),
          ]),
          results_by_covariate: { age: 'merged.csv' },
        },
        { model: 'developmental', status: 'complete', design: '~ age + juvenile_age', n_samples: 8,
          results_by_covariate: { juvenile_age: 'study.csv', age_25_vs_5: 'study.csv', adult_age: 'study.csv' } },
        ...[
          ['study_a', 'study.csv'],
          ['study_b', 'study-b.csv'],
          ['study_c', 'study-c.csv'],
        ].map(([study, file]) => ({
          model: 'per_study',
          study,
          status: 'complete',
          design: '~ age',
          n_samples: 8,
          results_by_covariate: { age: file },
        })),
      ],
    }),
  )
  await write(
    join(dir, 'merged.csv'),
    'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.12,0.001,0.0001,0.02,10,per year\nGENE2,-0.35,0.02,0.01,0.04,20,per year\nGENE3,0.23,0.5,0.4,0.05,7,per year\n',
  )
  await write(
    join(dir, 'reference_expression.csv'),
    'cell_type,gene,reference_bin,n_studies,n_samples,reference_cpm\nCD14 monocyte,GENE1,20-30,1,12,123.4\nCD14 monocyte,GENE2,20-30,1,12,56.7\n',
  )
  await write(
    join(dir, 'study.csv'),
    'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.44,0.004,0.001,0.03,9,per year\nGENE2,-0.21,0.03,0.02,0.05,19,per year\n',
  )
  await write(
    join(dir, 'study-b.csv'),
    'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.44,0.004,0.001,0.03,9,per year\nGENE2,-0.21,0.03,0.02,0.05,19,per year\n',
  )
  await write(
    join(dir, 'study-c.csv'),
    'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.44,0.004,0.001,0.03,9,per year\nGENE2,-0.21,0.03,0.02,0.05,19,per year\n',
  )
  await write(
    join(dir, 'fraction_model_diagnostics.tsv'),
    'cell_type\taifi_l1_parent\tresponse\tstudy\tn_samples\tmean_fraction\tresidual_sd\tbinomial_sampling_sd\tcovariates\nCD14 monocyte\tMonocyte\tfraction_within_aifi_l1_parent\tstudy_a\t12\t0.25\t0.12\t0.03\tage;sex\n',
  )
  await write(
    join(dir, 'fraction_model_coefficients.tsv'),
    'cell_type\taifi_l1_parent\tresponse\tstudy\tcovariate\tterm\tcontrast\testimate\tci_low\tci_high\tn_samples\nCD14 monocyte\tMonocyte\tfraction_within_aifi_l1_parent\tstudy_a\tage\tage\tper decade\t0.02\t0.01\t0.03\t12\n',
  )
  await write(
    join(dir, 'fraction_age_bins.tsv'),
    'cell_type\taifi_l1_parent\tstudy\tdenominator\tage_bin\treference_bin\tn_samples_bin\testimate\tci_low\tci_high\treference_mean_fraction\treference_adjusted_fraction\nCD14 monocyte\tMonocyte\tstudy_a\tl1_parent\t20-30\t20-30\t6\t0\t0\t0\t0.25\t0.25\nCD14 monocyte\tMonocyte\tstudy_a\tl1_parent\t30-40\t20-30\t6\t0.1\t0.05\t0.15\t\t\nCD14 monocyte\tMonocyte\tstudy_a\tretained_pbmc\t20-30\t20-30\t6\t0\t0\t0\t0.40\t0.40\nCD14 monocyte\tMonocyte\tstudy_a\tretained_pbmc\t30-40\t20-30\t6\t0.2\t0.15\t0.25\t\t\n',
  )
  await write(
    join(dir, 'report.html'),
    '<!doctype html><title>Fixture DE report</title><h1>Fixture DE report</h1>',
  )
  await write(join(trajectory, 'analysis_metadata.json'), JSON.stringify({ common_bins: bins }))
  const header = [
    'gene',
    'cell_type',
    'trajectory_id',
    'trajectory_cluster',
    'umap_1',
    'umap_2',
    'omnibus_fdr',
    ...bins,
    ...bins.map((b) => `raw_${b}`),
  ].join(',')
  await write(
    join(trajectory, 'cross_cell_type_trajectory_clusters.csv'),
    `${header}\nGENE1,CD14 monocyte,GENE1 | CD14 monocyte,1,1,2,0.001,0,1,2,3,4,5,6,0,0.1,0.2,0.3,0.4,0.5,0.6\nGENE2,CD14 monocyte,GENE2 | CD14 monocyte,2,2,3,0.02,0,-1,-2,-3,-4,-5,-6,0,-0.1,-0.2,-0.3,-0.4,-0.5,-0.7\nGENE2,Naive CD4 T cell,GENE2 | Naive CD4 T cell,2,3,4,0.09,0,1,2,3,4,5,6,0,0.1,0.2,0.3,0.4,0.5,0.6\n`,
  )
  await write(
    join(trajectory, 'report.html'),
    '<!doctype html><title>Fixture trajectory report</title><h1>Fixture trajectory report</h1>',
  )
  for (const [slug, label, hasFractions] of [
    ['naive-cd4-t-cell', 'Naive CD4 T cell', true],
    ['cd16-monocyte', 'CD16 monocyte', false],
  ]) {
    if (slug === 'naive-cd4-t-cell') await writeClock(slug, label, 'Study C')
    const other = join(input, slug, 'differential_expression')
    await write(
      join(other, 'run_metadata.json'),
      JSON.stringify({
        cell_type: label,
        analysis_mode: 'test_only',
        synthetic_test_data: true,
        interpretation_warning: 'Synthetic fixture: not biological evidence.',
        models: [
          {
            model: 'merged',
            status: 'complete',
            design: '~ age',
            n_samples: 10,
            results_by_covariate: { age: 'merged.csv' },
          },
          ...(hasFractions
            ? [
                {
                  model: 'developmental',
                  status: 'complete',
                  design: '~ age',
                  n_samples: 8,
                  results_by_covariate: { juvenile_age: 'study.csv' },
                },
              ]
            : []),
        ],
      }),
    )
    if (hasFractions)
      await write(
        join(other, 'study.csv'),
        'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.61,0.004,0.001,0.03,9,per decade\nGENE2,-0.67,0.03,0.02,0.05,19,per decade\n',
      )
    await write(
      join(other, 'reference_expression.csv'),
      `cell_type,gene,reference_bin,n_studies,n_samples,reference_cpm\n${label},GENE1,20-30,1,10,75\n${label},GENE2,20-30,1,10,25\n`,
    )
    await write(
      join(other, 'merged.csv'),
      `gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.11,0.01,0.001,0.02,10,per year\nGENE2,${hasFractions ? '-0.22,0.02' : '-0.32,0.03'},0.01,0.04,20,per year\n`,
    )
    if (hasFractions)
      await write(
        join(other, 'fraction_model_diagnostics.tsv'),
        `cell_type\taifi_l1_parent\tresponse\tstudy\tn_samples\tmean_fraction\tresidual_sd\tbinomial_sampling_sd\tcovariates\n${label}\tT cell\tfraction_within_aifi_l1_parent\tstudy_a\t10\t0.15\t0.1\t0.03\tage\n`,
      )
  }
  const studiesConfigPath = join(root, 'studies.json')
  await write(
    studiesConfigPath,
    JSON.stringify({
      studies: {
        study_a: {
          publication_doi: '10.1234/fixture',
          provenance: { celltypes: 'Fixture annotation', technology: 'Fixture technology' },
        },
      },
    }),
  )
  const publicDir = join(root, 'public')
  const manifest = await createRelease(input, join(publicDir, 'data'), {
    pathwayDirectory,
    pantherPath: join(root, 'no-panther'),
    studiesConfigPath,
  })
  return { publicDir, manifest }
}
