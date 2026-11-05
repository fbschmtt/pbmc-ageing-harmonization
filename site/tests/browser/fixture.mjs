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
  const analysis = join(input, 'cd14-monocyte/analysis')
  await write(join(analysis, 'report.json'), JSON.stringify({ status: 'complete', cell_type: 'CD14 monocyte', n_cells: 42, n_clusters: 2, embedding: { integration: { method: 'harmony' } } }))
  await write(join(analysis, 'markers.tsv'), 'names\tscores\tlogfoldchanges\tpvals_adj\tgroup\nTESTMARKER\t3.2\t2.4\t0.01\t0\nOTHER\t-2\t-1.5\t0.1\t0\nSECOND\t3\t1\t0.01\t1\n')
  await write(join(analysis, 'executed.ipynb'), JSON.stringify({ cells: [
    { cell_type: 'markdown', source: ['### Harmony-Integrated UMAP'] },
    { cell_type: 'code', outputs: [{ data: { 'image/png': 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/a1sAAAAASUVORK5CYII=' } }] },
    { cell_type: 'code', source: ['# cluster-fraction zeros included\ncluster_fractions = sample_cluster_fractions('], outputs: [
      { data: { 'text/markdown': ['_3 sample(s) excluded from cluster fractions because the split cell-type denominator has fewer than 10 cells._'] } },
      { data: { 'text/html': ['<details><summary>Show cluster-composition curves</summary><img src="data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/a1sAAAAASUVORK5CYII="></details>'] } },
      { data: { 'application/vnd.pbmc.cluster-fractions+json': { clusters: [
        { cluster: '0', image: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/a1sAAAAASUVORK5CYII=', fits: [{ study: null, samples: 3, meanFraction: 0.25, slopePpPerDecade: 1.5 }, { study: 'study_a', samples: 1, meanFraction: 0.5, slopePpPerDecade: null }] },
        { cluster: '1', image: 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/a1sAAAAASUVORK5CYII=', fits: [{ study: null, samples: 3, meanFraction: 0.75, slopePpPerDecade: -1.5 }] },
      ] } } },
    ] },
    { cell_type: 'code', source: ['# study-specific linear trends'], outputs: [] },
  ] }))
  const dir = join(input, 'cd14-monocyte/differential_expression')
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
          study_sample_counts: { study_a: 12 },
          results_by_covariate: { age: 'merged.csv' },
        },
        {
          model: 'per_study',
          study: 'study_a',
          status: 'complete',
          design: '~ age',
          n_samples: 8,
          results_by_covariate: { age: 'study.csv' },
        },
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
    'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.44,0.004,0.001,0.03,9,per decade\nGENE2,-0.21,0.03,0.02,0.05,19,per decade\n',
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
                  model: 'per_study',
                  study: 'study_a',
                  status: 'complete',
                  design: '~ age',
                  n_samples: 8,
                  results_by_covariate: { age: 'study.csv' },
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
      join(other, 'merged.csv'),
      'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.11,0.01,0.001,0.02,10,per year\nGENE2,-0.22,0.02,0.01,0.04,20,per year\n',
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
