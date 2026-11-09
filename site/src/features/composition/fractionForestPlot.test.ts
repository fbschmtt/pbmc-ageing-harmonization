import { expect, it } from 'vitest'
import type { FractionCoefficient } from '../../data/contracts.mjs'
import { clipForestInterval, fractionForestAxis } from './fractionForestPlot'

const coefficient = (overrides: Partial<FractionCoefficient>): FractionCoefficient => ({
  cellType: 'T cell',
  parent: 'T cell',
  study: 'study-a',
  covariate: 'age',
  term: 'age',
  contrast: 'per decade',
  estimate: 0.02,
  low: -0.2,
  high: 0.2,
  samples: 40,
  ...overrides,
})

it('sets the forest range from estimates and residual variation, not confidence interval tails', () => {
  const axis = fractionForestAxis(
    [coefficient({ low: -8, high: 8 })],
    [{ study: 'study-a', residualSd: 0.06, binomialSamplingSd: 0.01 }],
  )
  expect(axis.low).toBeGreaterThan(-1)
  expect(axis.high).toBeLessThan(1)
  expect(clipForestInterval({ low: -8, high: 8 }, axis.low, axis.high)).toMatchObject({
    low: axis.low,
    high: axis.high,
    clippedLow: true,
    clippedHigh: true,
  })
})
