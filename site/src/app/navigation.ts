export type Page =
  | 'differential'
  | 'overview'
  | 'datasets'
  | 'modelling'
  | 'genes'
  | 'trajectories'
  | 'geneTypeMatrix'
  | 'composition'
  | 'singleCell'
  | 'pathways'
  | 'about'
  | 'clock'

export type NavPage = { id: Page; label: string; description?: string }
export type NavigationItem = NavPage | {
  id: string
  label: string
  pages: NavPage[]
}

export const navigation: NavPage[] = [
  { id: 'overview', label: 'Overview' },
  { id: 'composition', label: 'Composition' },
  { id: 'differential', label: 'By cell type' },
  { id: 'genes', label: 'By gene' },
  { id: 'singleCell', label: 'Cell states' },
  { id: 'pathways', label: 'Pathways' },
  // Retain the former matrix hash as a route alias, but do not show it as a
  // separate navigation destination now that it is part of Pathways.
  { id: 'geneTypeMatrix', label: 'Gene × cell type' },
  { id: 'datasets', label: 'Datasets' },
  { id: 'modelling', label: 'Modelling' },
  { id: 'trajectories', label: 'Age-bin trajectories' },
  { id: 'clock', label: 'Transcriptomic clocks' },
  { id: 'about', label: 'About' },
]

const pageItem = (id: Page) => navigation.find((item) => item.id === id)!

export const navigationItems: NavigationItem[] = [
  pageItem('overview'),
  pageItem('composition'),
  {
    id: 'differential-menu',
    label: 'Differential expression',
    pages: [
      { id: 'differential', label: 'By cell type', description: 'Compare genes within a cell type' },
      { id: 'genes', label: 'By gene', description: 'Compare one gene across cell types' },
    ],
  },
  pageItem('singleCell'),
  {
    id: 'age-analysis-menu',
    label: 'More age analysis',
    pages: [
      { id: 'trajectories', label: 'Age-bin trajectories', description: 'Explore age-related gene profiles' },
      { id: 'clock', label: 'Transcriptomic clocks', description: 'Compare expression-based age predictions' },
    ],
  },
  pageItem('pathways'),
  pageItem('datasets'),
  pageItem('modelling'),
  pageItem('about'),
]

export function isPage(value: string): value is Page {
  return navigation.some((item) => item.id === value)
}
