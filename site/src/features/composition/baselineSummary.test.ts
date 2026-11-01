import { describe, expect, it } from 'vitest'
import { baselineSummary } from './baselineSummary'
import type { FractionAgeEstimate } from './FractionAgePlot'

function row(group: string, fraction: number | null, samples: number): FractionAgeEstimate {
  return {
    cellType: 'A',
    denominator: 'l1_parent',
    referenceBin: '20-30',
    ageBin: '20-30',
    group,
    estimate: 0,
    low: 0,
    high: 0,
    samples,
    referenceMeanFraction: fraction,
    referenceAdjustedFraction: fraction,
  }
}

describe('baseline summary', () => {
  it('weights each study equally rather than weighting by sample counts or repeated rows', () => {
    const a = row('a', 0.1, 5)
    expect(baselineSummary([a, a, row('b', 0.5, 50)])).toEqual({
      mean: 0.3,
      studies: [
        { study: 'a', value: 0.1, samples: 5 },
        { study: 'b', value: 0.5, samples: 50 },
      ],
    })
  })
  it('retains model estimates beyond physical bounds rather than clipping them', () => {
    expect(baselineSummary([row('a', 1.2, 5)])?.mean).toBe(1.2)
  })
  it('ignores missing baseline fractions and non-reference bins', () => {
    expect(baselineSummary([row('a', null, 5), { ...row('b', 0.4, 10), ageBin: '30-40' }])).toBeNull()
  })
})
