import { describe, expect, it } from 'vitest'
import { studySupportColor } from './studySupportScale'

describe('study support color scale', () => {
  it('uses viridis endpoints for zero and maximum study support', () => {
    expect(studySupportColor(0, 10)).toBe('rgb(68 1 84)')
    expect(studySupportColor(10, 10)).toBe('rgb(253 231 37)')
  })

  it('keeps every study count distinguishable on a ten-study scale', () => {
    expect(new Set(Array.from({ length: 11 }, (_, count) => studySupportColor(count, 10))).size).toBe(11)
  })

  it('clamps counts to the configured scale', () => {
    expect(studySupportColor(12, 10)).toBe(studySupportColor(10, 10))
  })
})
