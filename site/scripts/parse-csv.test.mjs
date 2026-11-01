import { describe, expect, it } from 'vitest'
import { parseCsv } from './parse-csv.mjs'

describe('CSV input', () => {
  it('preserves quoted commas, quotes and newlines', () => {
    expect(parseCsv('gene,note\r\nA,"line one,\nline ""two"""\r\n')).toEqual({
      headers: ['gene', 'note'],
      rows: [['A', 'line one,\nline "two"']],
    })
  })
  it('rejects malformed rows and incomplete quoted fields', () => {
    expect(() => parseCsv('gene,value\nA\n')).toThrow('wrong number')
    expect(() => parseCsv('gene,value\nA,"open')).toThrow('Unterminated')
  })
})
