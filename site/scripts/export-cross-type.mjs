import { readFile, writeFile, mkdir, rm } from 'node:fs/promises'
import { join, dirname } from 'node:path'
import { parseCsv } from './parse-csv.mjs'
import { parseTrajectoryMetadata } from './trajectory-metadata.mjs'
import { geneBucket } from '../src/data/geneBucket.mjs'

export const crossTypeSources = [
  'gene_recurrence.csv',
  'gene_cell_type_significance.csv',
  'gene_pattern_concordance.csv',
  'cross_cell_type_cluster_means.csv',
]
const number = (value) => (value == null || value.trim() === '' ? null : Number(value))
const boolean = (value) => (value === 'True' ? true : value === 'False' ? false : null)
async function table(directory, name, required) {
  const { headers, rows } = parseCsv(await readFile(join(directory, name), 'utf8'))
  if (required.some((key) => !headers.includes(key))) throw new Error(`${name}: missing cross-type columns`)
  return rows.map((values) => Object.fromEntries(headers.map((key, i) => [key, values[i]])))
}

export async function exportCrossType(input, output, cellTypes) {
  const directory = join(input, '..', 'combined', 'trajectory')
  const target = dirname(output)
  await rm(join(target, 'cross-type-index.json'), { force: true })
  await rm(join(target, 'cross-type-genes'), { recursive: true, force: true })
  try {
    await readFile(join(directory, 'gene_recurrence.csv'), 'utf8')
  } catch (error) {
    if (error.code === 'ENOENT') return
    throw error
  }
  const metadata = parseTrajectoryMetadata(await readFile(join(directory, 'analysis_metadata.json'), 'utf8'))
  const known = new Set(cellTypes.map((cell) => cell.label))
  const genes = await table(directory, 'gene_recurrence.csv', [
    'gene',
    'n_cell_types_tested',
    'n_cell_types_DE_significant',
  ])
  if (!genes.length) return
  const records = await table(directory, 'gene_cell_type_significance.csv', [
    'gene',
    'cell_type',
    'omnibus_fdr',
    'tested',
    'de_significant',
    'trajectory_cluster_significant',
    'largest_absolute_log2_fold_change',
    'largest_effect_age_bin',
  ])
  const patterns = await table(directory, 'gene_pattern_concordance.csv', [
    'gene',
    'n_cell_types_cluster_significant',
    'dominant_cluster_fraction',
  ])
  const bins = metadata.common_bins
  if (!Array.isArray(bins)) throw new Error('Cross-type metadata lacks common_bins')
  const clusters = await table(directory, 'cross_cell_type_cluster_means.csv', [
    'trajectory_cluster',
    'n_trajectories',
    ...bins,
  ])
  const details = new Map()
  for (const row of genes) {
    if (!row.gene || details.has(row.gene)) throw new Error(`Duplicate or empty cross-type gene: ${row.gene}`)
    details.set(row.gene, { gene: row.gene, records: [], pattern: null })
  }
  for (const row of records) {
    const detail = details.get(row.gene)
    if (
      !detail ||
      !known.has(row.cell_type) ||
      detail.records.some((record) => record.cellType === row.cell_type)
    )
      throw new Error(`Unknown or duplicate cross-type gene/cell type: ${row.gene}/${row.cell_type}`)
    detail.records.push({
      cellType: row.cell_type,
      padj: number(row.omnibus_fdr),
      tested: boolean(row.tested),
      significant: boolean(row.de_significant),
      clusterSignificant: boolean(row.trajectory_cluster_significant),
      effect: number(row.largest_absolute_log2_fold_change),
      ageBin: row.largest_effect_age_bin || null,
    })
  }
  for (const row of patterns) {
    const detail = details.get(row.gene)
    if (!detail || detail.pattern) throw new Error(`Unknown or duplicate cross-type pattern: ${row.gene}`)
    detail.pattern = {
      types: number(row.n_cell_types_cluster_significant),
      dominantFraction: number(row.dominant_cluster_fraction),
      meanDistance: number(row.mean_pairwise_shape_distance),
      sameClusterFraction: number(row.same_cluster_pair_fraction),
    }
  }
  const summaries = genes
    .map((row) => {
      const tested = number(row.n_cell_types_tested),
        significant = number(row.n_cell_types_DE_significant)
      const data = details.get(row.gene).records
      if (
        tested !== data.filter((record) => record.tested).length ||
        significant !== data.filter((record) => record.significant).length
      )
        throw new Error(`Cross-type recurrence disagrees with gene/type results: ${row.gene}`)
      return { gene: row.gene, tested, significant }
    })
    .sort((a, b) => b.significant - a.significant || b.tested - a.tested || a.gene.localeCompare(b.gene))
  const buckets = new Map(),
    geneBuckets = {}
  for (const [gene, detail] of details) {
    const bucket = geneBucket(gene)
    if (!buckets.has(bucket)) buckets.set(bucket, {})
    buckets.get(bucket)[gene] = detail
  }
  await mkdir(join(target, 'cross-type-genes'), { recursive: true })
  for (const [bucket, data] of buckets) {
    const path = `cross-type-genes/${bucket}.json`
    geneBuckets[bucket] = path
    await writeFile(join(target, path), JSON.stringify({ schemaVersion: 1, genes: data }))
  }
  const types = [...new Set(records.map((row) => row.cell_type))].sort()
  const exclusions = Object.entries(metadata.cell_types_excluded_from_cross_cell_type_analysis || {}).map(
    ([cellType, reason]) => ({ cellType, reason }),
  )
  const index = {
    schemaVersion: 1,
    cellTypes: types,
    genes: summaries,
    geneBuckets,
    deFdr: metadata.de_fdr_threshold,
    clusterFdr: metadata.cluster_fdr_threshold,
    commonBins: bins,
    referenceBin: `${metadata.reference_bin_start_age}-${metadata.reference_bin_start_age + metadata.bin_width_years}`,
    exclusions,
    clusters: clusters.map((row) => ({
      id: row.trajectory_cluster,
      count: number(row.n_trajectories),
      mean: bins.map((bin) => number(row[bin])),
    })),
  }
  // Use the same contracts at export, release verification, and browser loading.
  const { validAsset } = await import('../src/data/contracts.mjs')
  if (!validAsset('crossType', index)) throw new Error('Invalid cross-type summary asset')
  for (const data of buckets.values())
    if (!validAsset('crossTypeGenes', { schemaVersion: 1, genes: data }))
      throw new Error('Invalid cross-type gene asset')
  await writeFile(join(target, 'cross-type-index.json'), JSON.stringify(index))
}
