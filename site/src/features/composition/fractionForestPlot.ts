import type { FractionCoefficient } from '../../data/contracts.mjs'
import { rangeAxis } from '../../theme/axisScale'

export type FractionDiagnostic = {
  study: string
  residualSd: number | null
  binomialSamplingSd: number | null
}

export function fractionForestAxis(
  coefficients: FractionCoefficient[],
  diagnostics: FractionDiagnostic[],
) {
  const values = [
    0,
    ...coefficients.map((item) => item.estimate),
    ...diagnostics.flatMap((item) => [item.residualSd ?? 0, item.binomialSamplingSd ?? 0]),
  ].filter(Number.isFinite)
  const minimum = Math.min(...values)
  const maximum = Math.max(...values)
  const padding = Math.max((maximum - minimum) * 0.08, Math.max(Math.abs(minimum), Math.abs(maximum), 0.01) * 0.02)
  return rangeAxis(minimum, maximum, padding)
}

export function clipForestInterval(
  interval: { low: number; high: number },
  low: number,
  high: number,
) {
  return {
    low: Math.min(high, Math.max(interval.low, low)),
    high: Math.max(low, Math.min(interval.high, high)),
    clippedLow: interval.low < low,
    clippedHigh: interval.high > high,
  }
}
