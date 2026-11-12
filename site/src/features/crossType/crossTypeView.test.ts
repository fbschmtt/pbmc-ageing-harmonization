import { expect, it } from 'vitest'
import { validAsset, type CrossTypeIndex, type CrossTypeGene } from '../../data/contracts.mjs'
import { orderedResults, visibleGenes, patternExtent, patternPoints } from './crossTypeView'

const index: CrossTypeIndex = {
  schemaVersion: 1,
  cellTypes: ['A', 'B'],
  genes: [
    { gene: 'GENE1', significant: 2, tested: 2 },
    { gene: 'GENE2', significant: 0, tested: 1 },
  ],
  geneBuckets: {},
  deFdr: 0.05,
  clusterFdr: 0.001,
  referenceBin: '20-30',
  commonBins: [],
  exclusions: [],
  clusters: [],
}
it('filters by recurrence and gene name without changing tested denominators', () => {
  expect(visibleGenes(index, 'gene', 1)).toEqual([index.genes[0]])
  expect(visibleGenes(index, 'gene2', 0)[0].tested).toBe(1)
})
it('represents an absent cell type separately from a nonsignificant result', () => {
  const detail: CrossTypeGene = {
    gene: 'GENE1',
    pattern: null,
    records: [
      {
        cellType: 'A',
        tested: true,
        significant: false,
        clusterSignificant: false,
        padj: 0.5,
        effect: 0.2,
        ageBin: '30-40',
      },
    ],
  }
  expect(orderedResults(index, detail).map((row) => row.result?.tested ?? null)).toEqual([true, null])
})
it('rejects impossible significance and duplicate cell-type records', () => {
  const record = {
    cellType: 'A',
    tested: false,
    significant: true,
    clusterSignificant: false,
    padj: null,
    effect: null,
    ageBin: null,
  }
  expect(
    validAsset('crossTypeGenes', {
      schemaVersion: 1,
      genes: { GENE1: { gene: 'GENE1', records: [record], pattern: null } },
    }),
  ).toBe(false)
  record.significant = false
  expect(
    validAsset('crossTypeGenes', {
      schemaVersion: 1,
      genes: { GENE1: { gene: 'GENE1', records: [record, record], pattern: null } },
    }),
  ).toBe(false)
})
it('renders finite paths even for a zero profile or a single shared bin', () => {
  expect(patternPoints([0])).toBe('30,78')
  expect(patternPoints([0, -1, 1])).not.toContain('NaN')
})

it('rejects a malformed reference age bin', () => {
  expect(validAsset('crossType', { ...index, referenceBin: 'undefined-NaN' })).toBe(false)
})

it('uses one vertical scale without amplifying lower-amplitude pattern means', () => {
  const extent = patternExtent([{ mean: [0, 2] }, { mean: [0, 1] }])
  expect(extent).toBe(2)
  expect(patternPoints([0, 2], extent)).toBe('30,78 290,28')
  expect(patternPoints([0, 1], extent)).toBe('30,78 290,53')
  expect(patternExtent([])).toBe(0.25)
})
