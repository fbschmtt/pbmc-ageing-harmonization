import { describe, expect, it } from 'vitest'
import { formatAxisValue, positiveAxis, symmetricAxis } from './axisScale'

describe('axis scaling', () => {
  it('pads the domain without rounding it up to an outer tick', () => {
    const axis = symmetricAxis([-2.6, 1.4])
    expect(axis.extent).toBeCloseTo(2.73)
    expect(axis.ticks).toEqual([-2, -1, 0, 1, 2])
    expect(axis.extent).toBeGreaterThan(2.6)
    expect(Math.max(...axis.ticks)).toBeLessThan(2.6)
  })

  it('keeps fractional domains and ticks at their natural scale', () => {
    const axis = symmetricAxis([-0.43, 0.2])
    expect(axis.extent).toBeCloseTo(0.4515)
    expect(axis.ticks).toEqual([-0.4, -0.2, 0, 0.2, 0.4])
  })

  it('creates a positive axis whose data may extend above its last tick', () => {
    const axis = positiveAxis([0.12, 0.43])
    expect(axis.maximum).toBeCloseTo(0.4515)
    expect(axis.ticks).toEqual([0, 0.1, 0.2, 0.3, 0.4])
  })

  it('provides a stable fallback for zero-only data', () => {
    expect(symmetricAxis([0])).toEqual({ extent: 1, ticks: [-1, -0.5, 0, 0.5, 1] })
    expect(positiveAxis([0])).toEqual({ maximum: 1, ticks: [0, 0.2, 0.4, 0.6, 0.8, 1] })
  })

  it('formats whole, fractional, and very small axes without negative zero', () => {
    expect(formatAxisValue(10, 10)).toBe('10')
    expect(formatAxisValue(0.4, 0.45)).toBe('0.4')
    expect(formatAxisValue(0.006, 0.01)).toBe('0.006')
    expect(formatAxisValue(-Number.EPSILON, 0.5)).toBe('0')
  })
})
