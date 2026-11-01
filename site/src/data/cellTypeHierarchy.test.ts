import { describe, expect, it } from 'vitest'
import { compositionScopes, groupCellTypes } from './cellTypeHierarchy'

describe('groupCellTypes', () => {
  it('uses the declared biological hierarchy rather than input order', () => {
    const groups = groupCellTypes(['Memory B cell', 'Astro', 'EN_L2_3_IT'])
    expect(groups.map((group) => group.label)).toEqual([
      'Brain — neurons — excitatory',
      'Brain — glia',
      'PBMC — B and plasma cells',
    ])
  })

  it('retains unknown cell types in Other', () => {
    expect(groupCellTypes(['Future cell type'])).toEqual([
      { label: 'Other', cellTypes: ['Future cell type'] },
    ])
  })

  it('uses the same leaf definitions as the composition hierarchy', () => {
    expect(groupCellTypes(['CD8aa', 'PVM', 'SMC']).map((group) => group.label)).toEqual([
      'Brain — glia',
      'Brain — vascular and support',
      'PBMC — T cells',
    ])
  })
})

describe('compositionScopes', () => {
  it('provides nested lymphocyte scopes without mixing in myeloid cells', () => {
    const scopes = compositionScopes([
      { tissue: 'pbmc', cellType: 'Memory CD4 T cell' },
      { tissue: 'pbmc', cellType: 'Memory B cell' },
      { tissue: 'pbmc', cellType: 'CD56dim NK cell' },
      { tissue: 'pbmc', cellType: 'CD14 monocyte' },
    ])
    expect(scopes.find((scope) => scope.id === 'pbmc-lymphocytes')?.cellTypes).toEqual([
      'Memory CD4 T cell',
      'Memory B cell',
      'CD56dim NK cell',
    ])
    expect(scopes.find((scope) => scope.id === 'pbmc-myeloid')?.cellTypes).toEqual(['CD14 monocyte'])
  })

  it('keeps unclassified types visible at tissue level but out of narrower scopes', () => {
    const scopes = compositionScopes([
      { tissue: 'pbmc', cellType: 'Memory B cell' },
      { tissue: 'pbmc', cellType: 'Future cell type' },
    ])
    expect(scopes.find((scope) => scope.id === 'pbmc')?.cellTypes).toEqual([
      'Memory B cell',
      'Future cell type',
    ])
    expect(scopes.find((scope) => scope.id === 'pbmc-lymphocytes')?.cellTypes).toEqual(['Memory B cell'])
  })
})
