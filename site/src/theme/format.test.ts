import { describe, expect, it } from 'vitest'
import { formatCpm, formatNumber } from './format'

describe('shared number formatting', () => {
  it('uses two decimal places for values with magnitude at least one', () => {
    expect(formatNumber(12)).toBe('12.00')
    expect(formatNumber(-1.237)).toBe('-1.24')
  })

  it('uses two significant digits below one and keeps zero readable', () => {
    expect(formatNumber(0.456)).toBe('0.46')
    expect(formatNumber(0.00456)).toBe('0.0046')
    expect(formatNumber(6.89e-47)).toBe('6.9e-47')
    expect(formatNumber(0)).toBe('0')
  })

  it('handles missing values and applies the rule to CPM labels', () => {
    expect(formatNumber(null)).toBe('—')
    expect(formatNumber(Number.NaN)).toBe('—')
    expect(formatCpm(0.00456)).toBe('0.0046 CPM')
  })
})
