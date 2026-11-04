import { expect, it } from 'vitest'
import type { FractionAgeEstimate } from '../../data/contracts.mjs'
import { absoluteFractionPoints, clipFractionInterval, pointBasedFractionAxis } from './ageFractionPlot'

const row = (overrides: Partial<FractionAgeEstimate>): FractionAgeEstimate => ({
  cellType: 'Treg',
  denominator: 'l1_parent',
  referenceBin: '20-30',
  referenceMeanFraction: null,
  referenceAdjustedFraction: null,
  ageBin: '20-30',
  group: 'study-a',
  estimate: 0,
  low: 0,
  high: 0,
  samples: 8,
  ...overrides,
})

it('anchors each study trajectory at its adjusted absolute reference fraction', () => {
  const points = absoluteFractionPoints([
    row({ referenceAdjustedFraction: 0.25 }),
    row({ ageBin: '30-40', estimate: 0.1, low: 0.05, high: 0.15 }),
    row({ group: 'study-b', referenceAdjustedFraction: 0.4 }),
    row({ group: 'study-b', ageBin: '30-40', estimate: -0.1, low: -0.2, high: 0 }),
  ])
  expect(points).toHaveLength(4)
  expect(points[0].estimate).toBeCloseTo(0.25)
  expect(points[1].estimate).toBeCloseTo(0.35)
  expect(points[1].low).toBeCloseTo(0.3)
  expect(points[1].high).toBeCloseTo(0.4)
  expect(points[2].estimate).toBeCloseTo(0.4)
  expect(points[3].estimate).toBeCloseTo(0.3)
  expect(points[3].low).toBeCloseTo(0.2)
  expect(points[3].high).toBeCloseTo(0.4)
})

it('sets the axis from estimates alone and marks confidence intervals past its limits', () => {
  const points = absoluteFractionPoints([
    row({ referenceAdjustedFraction: 0.2, low: -20, high: 20 }),
    row({ ageBin: '30-40', estimate: 0.2, low: -1, high: 1 }),
  ])
  const axis = pointBasedFractionAxis(points)
  expect(axis.low).toBeCloseTo(0.18)
  expect(axis.high).toBeCloseTo(0.42)
  expect(clipFractionInterval(points[1], axis.low, axis.high)).toMatchObject({
    low: axis.low,
    high: axis.high,
    clippedLow: true,
    clippedHigh: true,
  })
})

it('omits trajectories without an adjusted reference estimate', () => {
  expect(absoluteFractionPoints([row({ referenceAdjustedFraction: null })])).toEqual([])
})
