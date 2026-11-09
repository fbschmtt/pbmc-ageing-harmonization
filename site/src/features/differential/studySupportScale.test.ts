import { describe, expect, it } from 'vitest'
import { studySupportColor } from './studySupportScale'

describe('study support color scale', () => {
  it('keeps zero support neutral and maps support onto viridis', () => {
    expect(studySupportColor(0, 10)).toBe('#b9bdbc')
    expect(studySupportColor(1, 10)).toBe('#440154')
    expect(studySupportColor(10, 10)).toBe('#fde725')
  })

  it('keeps every study count distinguishable on a ten-study scale', () => {
    expect(new Set(Array.from({ length: 11 }, (_, count) => studySupportColor(count, 10))).size).toBe(11)
  })

  it('clamps counts to the configured scale', () => {
    expect(studySupportColor(12, 10)).toBe(studySupportColor(10, 10))
  })
})
