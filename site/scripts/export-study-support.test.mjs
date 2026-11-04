import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { afterEach, expect, it } from 'vitest'
import { addStudySupport } from './export-study-support.mjs'

const roots = []
afterEach(async () => { await Promise.all(roots.splice(0).map((root) => rm(root, { recursive: true, force: true }))) })

it('counts finite study p-values for each matching contrast, including nonsignificant lookup chunks', async () => {
  const root = await mkdtemp(join(tmpdir(), 'study-support-')); roots.push(root)
  await mkdir(join(root, 'de-models'))
  const models = [
    { id: 'merged', model: 'merged', study: null, covariate: 'age', contrast: 'per year', recordsAsset: 'de-models/merged.json', recordBuckets: { aa: 'de-models/lookup.json' } },
    { id: 'a', model: 'per_study', study: 'a', covariate: 'age', contrast: 'per year', recordsAsset: 'de-models/a.json', recordBuckets: { aa: 'de-models/a-other.json' } },
    { id: 'b', model: 'per_study', study: 'b', covariate: 'age', contrast: 'per year', recordsAsset: 'de-models/b.json' },
    { id: 'different', model: 'per_study', study: 'c', covariate: 'age', contrast: 'different', recordsAsset: 'de-models/c.json' },
  ]
  const records = {
    merged: [{ gene: 'GENE', padj: 0.01 }, { gene: 'MISSING', padj: 0.01 }],
    lookup: [{ gene: 'OTHER', padj: 0.5 }],
    a: [{ gene: 'GENE', padj: 0.01 }],
    'a-other': [{ gene: 'OTHER', padj: 0.05 }, { gene: 'MISSING', padj: null }],
    b: [{ gene: 'GENE', padj: 0.04 }, { gene: 'OTHER', padj: 0.001 }],
    c: [{ gene: 'GENE', padj: 0.001 }],
  }
  await writeFile(join(root, 'de-index.json'), JSON.stringify({ cellTypes: [{ id: 'cell', models }] }))
  for (const [name, rows] of Object.entries(records))
    await writeFile(join(root, `de-models/${name}.json`), JSON.stringify({ cellTypes: [{ models: [{ records: rows }] }] }))
  await addStudySupport(root)
  const read = async (name) => JSON.parse(await readFile(join(root, `de-models/${name}.json`), 'utf8')).cellTypes[0].models[0].records
  expect(await read('merged')).toMatchObject([
    { gene: 'GENE', studySupport: { significant: 2, tested: 2 } },
    { gene: 'MISSING', studySupport: null },
  ])
  expect(await read('lookup')).toMatchObject([{ gene: 'OTHER', studySupport: { significant: 1, tested: 2 } }])
})
