import { afterEach, describe, expect, it } from 'vitest'
import { mkdtemp, mkdir, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { compareAnalysisMetadata, releaseContext } from './release-provenance.mjs'

let root
afterEach(async () => {
  if (root) await rm(root, { recursive: true, force: true })
})
describe('release provenance', () => {
  it('compares model metadata and detects contradictions without claiming file lineage', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-provenance-'))
    const input = join(root, 'analysis/per_type'),
      local = join(input, 'mono/differential_expression'),
      combined = join(root, 'analysis/combined/differential_expression')
    await mkdir(local, { recursive: true })
    await mkdir(combined, { recursive: true })
    const model = {
      model: 'merged',
      status: 'complete',
      study: null,
      design: '~ age',
      n_samples: 10,
      covariates: ['age'],
      results_by_covariate: { age: 'merged.csv' },
    }
    await writeFile(
      join(local, 'run_metadata.json'),
      JSON.stringify({
        cell_type: 'Monocyte',
        analysis_mode: 'standard',
        sample_inclusion_applied: { a: 1, b: 2 },
        models: [model],
      }),
    )
    const manifest = {
      analysis_mode: 'standard',
      sample_inclusion_applied: { b: 2, a: 1 },
      cell_types: [
        {
          slug: 'mono',
          cell_type: 'Monocyte',
          models: [
            {
              ...model,
              results_by_covariate: { age: '../../per_type/mono/differential_expression/merged.csv' },
            },
          ],
        },
      ],
    }
    const path = join(combined, 'differential_expression.json')
    await writeFile(path, JSON.stringify(manifest))
    const result = await compareAnalysisMetadata(input, [{ id: 'mono' }])
    expect(result.status).toBe('consistent')
    expect(result.note).toContain('not established')
    manifest.cell_types[0].models[0].design = '~ sex'
    await writeFile(path, JSON.stringify(manifest))
    await expect(compareAnalysisMetadata(input, [{ id: 'mono' }])).rejects.toThrow('design disagrees')
  })
  it('reports only observed studies and labels missing references as unknown', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-context-'))
    const input = join(root, 'analysis/per_type')
    await mkdir(input, { recursive: true })
    const config = join(root, 'studies.json')
    await writeFile(
      config,
      JSON.stringify({
        studies: {
          one: { publication_doi: '10.1234/fixture', provenance: { celltypes: 'Fixture annotation' } },
          unused: { publication_doi: '10.1234/unused' },
        },
      }),
    )
    const context = await releaseContext(
      input,
      { cellTypes: [{ models: [{ study: 'one', studies: { unknown: 5 } }] }] },
      config,
    )
    expect(context.studies.map((s) => s.id)).toEqual(['one', 'unknown'])
    expect(context.studies[0].publicationDoi).toBe('10.1234/fixture')
    expect(context.studies[1].publicationDoi).toBeNull()
    expect(context.studyReference.note).toContain('not a verified')
  })
  it('includes harmonization studies independently from DE studies', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-mixed-context-'))
    const input = join(root, 'analysis/per_type')
    await mkdir(input, { recursive: true })
    const context = await releaseContext(
      input,
      { cellTypes: [{ models: [{ study: 'synthetic', studies: {} }] }] },
      join(root, 'missing.json'),
      [{ kind: 'report', study: 'real-study' }],
    )
    expect(context.studies).toEqual([
      expect.objectContaining({ id: 'real-study', inModels: false, harmonizationAvailable: true }),
      expect.objectContaining({ id: 'synthetic', inModels: true, harmonizationAvailable: false }),
    ])
  })
})
