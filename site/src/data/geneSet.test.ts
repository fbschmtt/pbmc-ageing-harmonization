import { describe, expect, it } from 'vitest'
import { parseGeneSet } from './geneSet'

describe('parseGeneSet', () => {
  it('accepts commas and all whitespace without regular expressions', () => {
    expect([...parseGeneSet(' lamp1,ctsd\nGBA\tNPC1 ')]).toEqual(['LAMP1', 'CTSD', 'GBA', 'NPC1'])
  })

  it('deduplicates symbols after normalizing case', () => {
    expect([...parseGeneSet('Lamp1 LAMP1 lamp1')]).toEqual(['LAMP1'])
  })
})
