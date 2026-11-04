export type Page =
  | 'differential'
  | 'overview'
  | 'datasets'
  | 'modelling'
  | 'genes'
  | 'cellTypes'
  | 'trajectories'
  | 'geneTypeMatrix'
  | 'composition'
  | 'singleCell'
  | 'pathways'
  | 'download'

export const navigation: { id: Page; label: string }[] = [
  { id: 'overview', label: 'Overview' },
  { id: 'differential', label: 'Differential expression' },
  { id: 'composition', label: 'Composition' },
  { id: 'genes', label: 'Genes' },
  { id: 'cellTypes', label: 'Cell types' },
  { id: 'trajectories', label: 'Cross-type trajectory map' },
  { id: 'geneTypeMatrix', label: 'Gene × cell type' },
  { id: 'singleCell', label: 'Single-cell' },
  { id: 'pathways', label: 'Pathways' },
  { id: 'datasets', label: 'Datasets' },
  { id: 'modelling', label: 'Modelling' },
  { id: 'download', label: 'Download' },
]
