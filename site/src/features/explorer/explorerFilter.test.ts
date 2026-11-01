import { describe, expect, it } from 'vitest'
import type { TrajectoryRecord } from '../../data/atlas'
import {
  compileFocusQuery,
  hasRegexSpecialCharacters,
  maxFocusPatternLength,
  recordMatchesFocus,
} from './explorerFilter'

const record = {
  gene: 'MBP',
  cellType: 'Oligo',
} as TrajectoryRecord

describe('recordMatchesFocus', () => {
  it('applies anchored regular expressions independently to cell types', () => {
    expect(recordMatchesFocus(record, '^Oligo', /^Oligo/i)).toBe(true)
  })

  it('applies anchored regular expressions independently to genes', () => {
    expect(recordMatchesFocus(record, '^MBP$', /^MBP$/i)).toBe(true)
  })

  it('detects input that will receive regular-expression semantics', () => {
    expect(hasRegexSpecialCharacters('Oligo')).toBe(false)
    expect(hasRegexSpecialCharacters('^Oligo')).toBe(true)
    expect(hasRegexSpecialCharacters('CD8+ T cell')).toBe(true)
  })
})

describe('compileFocusQuery', () => {
  it('compiles ordinary and anchored input case-insensitively', () => {
    expect(compileFocusQuery('oligo').expression?.test('Oligo')).toBe(true)
    expect(compileFocusQuery('^Oligo$').expression?.test('Oligo')).toBe(true)
  })

  it('reports malformed expressions without producing a matcher', () => {
    expect(compileFocusQuery('(')).toEqual({
      expression: null,
      error: 'This regular expression is not valid.',
    })
  })

  it('rejects patterns beyond the supported length', () => {
    const result = compileFocusQuery('A'.repeat(maxFocusPatternLength + 1))
    expect(result.expression).toBeNull()
    expect(result.error).toContain(`${maxFocusPatternLength} characters`)
  })
})
