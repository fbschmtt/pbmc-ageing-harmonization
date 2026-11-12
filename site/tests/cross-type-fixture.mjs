import { mkdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'

export async function writeCrossTypeFixture(directory) {
  await mkdir(directory, { recursive: true })
  const bins = ['20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90']
  const write = (name, text) => writeFile(join(directory, name), text)
  await write(
    'analysis_metadata.json',
    JSON.stringify({
      common_bins: bins,
      de_fdr_threshold: 0.05,
      cluster_fdr_threshold: 0.001,
      reference_bin_start_age: 20,
      bin_width_years: 10,
      cell_types_excluded_from_cross_cell_type_analysis: {
        'CD16 monocyte': 'Fewer than seven supported age bins',
      },
    }),
  )
  await write(
    'gene_recurrence.csv',
    'gene,n_cell_types_tested,n_cell_types_DE_significant,significant_cell_types\nGENE1,2,2,CD14 monocyte; Naive CD4 T cell\nGENE2,2,1,CD14 monocyte\nGENE3,1,0,\n',
  )
  await write(
    'gene_cell_type_significance.csv',
    'gene,cell_type,omnibus_fdr,tested,de_significant,trajectory_cluster_significant,largest_absolute_log2_fold_change,largest_effect_age_bin\nGENE1,CD14 monocyte,0.0001,True,True,True,0.7,70-80\nGENE1,Naive CD4 T cell,0.0002,True,True,True,-0.4,80-90\nGENE2,CD14 monocyte,0.03,True,True,False,-0.7,80-90\nGENE2,Naive CD4 T cell,0.2,True,False,False,0.2,40-50\nGENE3,CD14 monocyte,,False,False,False,,\nGENE3,Naive CD4 T cell,0.4,True,False,False,0.1,30-40\n',
  )
  await write(
    'gene_pattern_concordance.csv',
    'gene,n_cell_types_cluster_significant,dominant_cluster_fraction,mean_pairwise_shape_distance,same_cluster_pair_fraction\nGENE1,2,0.5,2.8,0\n',
  )
  await write(
    'cross_cell_type_cluster_means.csv',
    `trajectory_cluster,n_trajectories,${bins.join(',')}\n1,2,0,0.2,0.4,0.6,0.8,1,1.2\n2,3,0,-0.2,-0.4,-0.6,-0.8,-1,-1.2\n`,
  )
}
