import { describe, expect, it } from 'vitest'
import { validAsset } from './contracts.mjs'

describe('release contracts', () => {
  const gene = { gene: 'A', effect: 0.2, padj: 0.01, pvalue: null, se: null, baseMean: null }
  const model = {
    id: 'merged:all:age',
    model: 'merged',
    study: null,
    covariate: 'age',
    contrast: 'per year',
    design: '~ age',
    samples: 10,
    studies: { A: 10 },
    records: [gene],
  }
  const cell = {
    id: 'a',
    label: 'A',
    tissue: 'pbmc',
    analysisMode: 'standard',
    warning: null,
    trajectoryAsset: null,
    models: [model],
  }
  const release = { schemaVersion: 1, source: 'pseudobulk', cellTypes: [cell] }
  it('accepts valid data and null estimates', () => expect(validAsset('de', release)).toBe(true))
  it('rejects incompatible versions and malformed nested records', () => {
    expect(validAsset('de', { ...release, schemaVersion: 2 })).toBe(false)
    expect(
      validAsset('de', {
        ...release,
        cellTypes: [{ ...cell, models: [{ ...model, records: [{ ...gene, padj: 2 }] }] }],
      }),
    ).toBe(false)
    expect(
      validAsset('preview', {
        schemaVersion: 1,
        referenceBin: null,
        design: null,
        genes: { A: [{ bin: '20-30', effect: '0' }] },
      }),
    ).toBe(false)
  })
})
