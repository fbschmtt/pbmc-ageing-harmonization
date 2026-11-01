import { describe, expect, it } from 'vitest'
import type { Cell, Gene, Model } from '../../data/differential'
import { filterGenes, hasPlottableEstimate, initialCell, visibleGenes } from './differentialView'

const gene = (name: string, effect: number | null, padj: number | null): Gene => ({
  gene: name,
  effect,
  padj,
  pvalue: null,
  se: null,
  baseMean: null,
})
const model = (records: Gene[]): Model => ({
  id: 'merged:all:age',
  model: 'merged',
  study: null,
  covariate: 'age',
  contrast: 'per year',
  design: '~ age',
  samples: 10,
  studies: {},
  records,
})
const cell = (label: string, records: Gene[]): Cell => ({
  id: label,
  label,
  tissue: 'pbmc',
  analysisMode: 'standard',
  warning: null,
  trajectoryAsset: null,
  models: [model(records)],
})
describe('differential view selection', () => {
  it('shows significant records and exact requested genes for older unfiltered exports', () => {
    const records = [gene('SIG', 1, 0.01), gene('EDGE', 1, 0.05), gene('OTHER', 1, 0.7)]
    expect(visibleGenes(records, new Set()).map((record) => record.gene)).toEqual(['SIG'])
    expect(visibleGenes(records, new Set(['OTHER'])).map((record) => record.gene)).toEqual(['SIG', 'OTHER'])
  })
  it('uses the preferred cell before looking for the preferred gene', () => {
    const cells = [cell('A', [gene('X', 1, 0.1)]), cell('B', [gene('Y', 1, 0.01)])]
    expect(initialCell(cells, 'Y', 'A').label).toBe('A')
    expect(initialCell(cells, 'Y', null).label).toBe('B')
  })
  it('filters to finite estimates and sorts by adjusted p value', () => {
    const records = [gene('ZZ', 1, 0.2), gene('ZA', 1, 0.01), gene('ZX', null, 0.001), gene('BAD', 1, 2)]
    expect(filterGenes(model(records), ' z ').map((x) => x.gene)).toEqual(['ZA', 'ZZ'])
    expect(hasPlottableEstimate(records[2])).toBe(false)
    expect(hasPlottableEstimate(records[3])).toBe(false)
  })
})
