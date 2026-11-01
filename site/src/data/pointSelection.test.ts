import { describe, expect, it } from 'vitest'
import { choosePoint } from './pointSelection'

const points = [
  { gene: 'A', p: 0.1 },
  { gene: 'B', p: 0.001 },
  { gene: 'C', p: null },
]

describe('initial point selection', () => {
  it('prefers a saved plottable gene over a more significant result', () => {
    expect(
      choosePoint(
        points,
        (point) => point.p,
        (point) => point.gene === 'A',
      )?.gene,
    ).toBe('A')
  })
  it('uses the most significant point when the saved gene is unavailable', () => {
    expect(
      choosePoint(
        points,
        (point) => point.p,
        (point) => point.gene === 'missing',
      )?.gene,
    ).toBe('B')
  })
  it('does not choose missing, nonfinite or invalid p-values', () => {
    const invalid = [{ p: null }, { p: NaN }, { p: -1 }, { p: 2 }]
    expect(choosePoint(invalid, (point) => point.p)).toBeNull()
  })
  it('falls back when a preferred record has no usable p-value', () => {
    expect(
      choosePoint(
        points,
        (point) => point.p,
        (point) => point.gene === 'C',
      )?.gene,
    ).toBe('B')
  })
})
