import { expect, it } from 'vitest'
import type { FractionIndex } from '../../data/fractions'
import { fractionView } from './fractionView'

it('selects the requested denominator and its skipped statuses', () => {
  const common = {
    cellType: 'A',
    referenceBin: '20-30',
    referenceMeanFraction: null,
    ageBin: '20-30',
    group: 'study',
    estimate: 0,
    low: 0,
    high: 0,
    samples: 5,
  }
  const data: FractionIndex = {
    schemaVersion: 1,
    records: [
      {
        cellType: 'A',
        parent: 'P',
        study: 'study',
        samples: 5,
        meanFraction: 0.2,
        residualSd: null,
        binomialSamplingSd: null,
        covariates: '',
        warning: null,
      },
    ],
    coefficients: [],
    ageEstimates: [
      { ...common, denominator: 'l1_parent', referenceAdjustedFraction: 0.25 },
      { ...common, denominator: 'retained_pbmc', referenceAdjustedFraction: 0.4 },
    ],
    ageStatus: [{ cellType: 'A', denominator: 'l1_parent', status: 'skipped', reason: 'unsupported' }],
  }
  expect(fractionView(data, 'A', 'l1_parent').baseline?.mean).toBe(0.25)
  expect(fractionView(data, 'A', 'l1_parent').skipped).toHaveLength(1)
  expect(fractionView(data, 'A', 'retained_pbmc').baseline?.mean).toBe(0.4)
  expect(fractionView(data, 'A', 'retained_pbmc').skipped).toHaveLength(0)
})
