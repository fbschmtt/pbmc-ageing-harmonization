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
  it('validates reference-bin expression values and support counts', () => {
    expect(validAsset('referenceExpression', {
      schemaVersion: 1,
      cellType: 'Naive CD8 T cell',
      referenceBin: '20-30',
      method: 'equal-study mean CPM',
      genes: { GENE1: { cpm: 12.5, studies: 2, samples: 8 } },
    })).toBe(true)
    expect(validAsset('referenceExpression', {
      schemaVersion: 1,
      cellType: 'Naive CD8 T cell',
      referenceBin: '20-30',
      method: 'equal-study mean CPM',
      genes: { GENE1: { cpm: -1, studies: 2, samples: 8 } },
    })).toBe(false)
  })
  it('validates compact trajectory map points separately from profile records', () => {
    const map = {
      schemaVersion: 1,
      release: 'Current pipeline',
      records: [{ id: 'G1 | T cell', x: 0.2, y: -0.1, cellType: 'T cell', tissue: 'pbmc', gene: 'G1', adjustedPValue: 0.01, cluster: '1' }],
      tissues: ['pbmc'],
      cellTypes: ['T cell'],
      genes: ['G1'],
      geneBuckets: { g1: 'trajectory-details/g1.json' },
    }
    expect(validAsset('trajectoryMap', map)).toBe(true)
    expect(validAsset('trajectoryMap', { ...map, records: [{ ...map.records[0], x: 'invalid' }] })).toBe(false)
  })
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
