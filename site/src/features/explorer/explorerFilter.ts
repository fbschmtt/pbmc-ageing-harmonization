import type { TrajectoryPoint } from '../../data/atlas'

export const maxFocusPatternLength = 160

export function compileFocusQuery(query: string) {
  if (!query) return { expression: null, error: null }
  if (query.length > maxFocusPatternLength) {
    return {
      expression: null,
      error: `Regex patterns are limited to ${maxFocusPatternLength} characters.`,
    }
  }
  try {
    return { expression: new RegExp(query, 'i'), error: null }
  } catch {
    return { expression: null, error: 'This regular expression is not valid.' }
  }
}

export function recordMatchesFocus(record: TrajectoryPoint, query: string, expression: RegExp | null) {
  if (!query) return true
  return Boolean(expression?.test(record.gene) || expression?.test(record.cellType))
}

export function hasRegexSpecialCharacters(query: string) {
  return /[\\^$.*+?()[\]{}|]/.test(query)
}
