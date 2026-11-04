import type { FractionAgeEstimate } from '../../data/contracts.mjs'
import { rangeAxis } from '../../theme/axisScale'

export type AbsoluteFractionPoint = {
  item: FractionAgeEstimate
  estimate: number
  low: number
  high: number
}

export function absoluteFractionPoints(estimates: FractionAgeEstimate[]): AbsoluteFractionPoint[] {
  const referenceByGroup = new Map<string, number>()
  for (const item of estimates) {
    const value = item.referenceAdjustedFraction
    if (item.ageBin === item.referenceBin && value != null && Number.isFinite(value))
      referenceByGroup.set(item.group, value)
  }
  return estimates.flatMap((item) => {
    const reference = referenceByGroup.get(item.group)
    if (reference === undefined || ![item.estimate, item.low, item.high].every(Number.isFinite)) return []
    return [{
      item,
      estimate: reference + item.estimate,
      low: reference + item.low,
      high: reference + item.high,
    }]
  })
}

export function pointBasedFractionAxis(points: AbsoluteFractionPoint[]) {
  const values = points.map((item) => item.estimate)
  const minimum = Math.min(...values)
  const maximum = Math.max(...values)
  const padding = Math.max((maximum - minimum) * 0.1, 0.01)
  return rangeAxis(minimum, maximum, padding)
}

export function clipFractionInterval(point: AbsoluteFractionPoint, low: number, high: number) {
  return {
    low: Math.min(high, Math.max(point.low, low)),
    high: Math.max(low, Math.min(point.high, high)),
    clippedLow: point.low < low,
    clippedHigh: point.high > high,
  }
}
