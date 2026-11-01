import { describe, expect, it } from 'vitest'
import { selectInformativeTerms } from './gene-metadata-terms.mjs'

describe('selectInformativeTerms', () => {
  it('prefers rarer terms and excludes GO roots', () => {
    const terms = new Map([
      ['GO:0008150', { id: 'GO:0008150', label: 'biological process' }],
      ['GO:0005622', { id: 'GO:0005622', label: 'intracellular anatomical structure' }],
      ['GO:0006810', { id: 'GO:0006810', label: 'transport' }],
      ['GO:0006888', { id: 'GO:0006888', label: 'ER to Golgi vesicle-mediated transport' }],
    ])
    const frequencies = new Map([
      ['GO:0008150', 100],
      ['GO:0005622', 2],
      ['GO:0006810', 60],
      ['GO:0006888', 4],
    ])

    expect(selectInformativeTerms(terms, frequencies, 2)).toEqual([
      { id: 'GO:0006888', label: 'ER to Golgi vesicle-mediated transport' },
      { id: 'GO:0006810', label: 'transport' },
    ])
  })
})
