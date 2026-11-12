import { afterEach, expect, it } from 'vitest'
import { mkdtemp, mkdir, readFile, writeFile, rm } from 'node:fs/promises'
import { join } from 'node:path'
import { tmpdir } from 'node:os'
import { exportCrossType } from './export-cross-type.mjs'
import { sourceInventory } from './release-sources.mjs'
import { geneBucket } from '../src/data/geneBucket.mjs'
import { writeCrossTypeFixture } from '../tests/cross-type-fixture.mjs'

const roots = []
afterEach(async () => {
  for (const root of roots.splice(0)) await rm(root, { recursive: true, force: true })
})
async function fixture() {
  const root = await mkdtemp(join(tmpdir(), 'cross-type-'))
  roots.push(root)
  const input = join(root, 'analysis/per_type'),
    output = join(root, 'data/de-index.json')
  await mkdir(input, { recursive: true })
  await mkdir(join(root, 'data'))
  const combined = join(input, '../combined/trajectory')
  await writeCrossTypeFixture(combined)
  return { input, output, combined, root, cells: [{ label: 'CD14 monocyte' }, { label: 'Naive CD4 T cell' }] }
}
it('exports named aggregate tables with lazy gene details and preserves signed effects', async () => {
  const { input, output, root, cells } = await fixture()
  await exportCrossType(input, output, cells)
  const index = JSON.parse(await readFile(join(root, 'data/cross-type-index.json'), 'utf8'))
  expect(index.genes[0]).toEqual({ gene: 'GENE1', significant: 2, tested: 2 })
  const bucket = JSON.parse(
    await readFile(join(root, 'data', index.geneBuckets[geneBucket('GENE1')]), 'utf8'),
  )
  expect(bucket.genes.GENE1.records[1].effect).toBe(-0.4)
  expect(bucket.genes.GENE1.pattern.sameClusterFraction).toBe(0)
  expect(index.clusters[0].mean).toHaveLength(7)
  expect(index.referenceBin).toBe('20-30')
  expect(index.exclusions[0].cellType).toBe('CD16 monocyte')
  const sources = await sourceInventory(input)
  expect(sources.some((source) => source.path.endsWith('/gene_cell_type_significance.csv'))).toBe(true)
})
it('keeps untested results distinct from nonsignificant tests', async () => {
  const { input, output, root, cells } = await fixture()
  await exportCrossType(input, output, cells)
  const bucket = JSON.parse(
    await readFile(join(root, `data/cross-type-genes/${geneBucket('GENE3')}.json`), 'utf8'),
  )
  expect(bucket.genes.GENE3.records.map((row) => [row.tested, row.significant, row.padj])).toEqual([
    [false, false, null],
    [true, false, 0.4],
  ])
})
it('allows older outputs without the optional summary and removes stale optional assets', async () => {
  const { input, output, root, combined, cells } = await fixture()
  await exportCrossType(input, output, cells)
  await rm(join(combined, 'gene_recurrence.csv'))
  await exportCrossType(input, output, cells)
  await expect(readFile(join(root, 'data/cross-type-index.json'))).rejects.toThrow()
})
it('rejects partial named artifacts rather than silently fabricating results', async () => {
  const { input, output, combined, cells } = await fixture()
  await rm(join(combined, 'gene_cell_type_significance.csv'))
  await expect(exportCrossType(input, output, cells)).rejects.toThrow()
})
it('rejects recurrence counts that disagree with the underlying gene/type tests', async () => {
  const { input, output, combined, cells } = await fixture()
  const path = join(combined, 'gene_recurrence.csv')
  await writeFile(path, (await readFile(path, 'utf8')).replace('GENE1,2,2', 'GENE1,2,1'))
  await expect(exportCrossType(input, output, cells)).rejects.toThrow('disagrees')
})
it('rejects duplicate gene/type records', async () => {
  const { input, output, combined, cells } = await fixture()
  const path = join(combined, 'gene_cell_type_significance.csv')
  await writeFile(
    path,
    (await readFile(path, 'utf8')) + 'GENE1,CD14 monocyte,0.0001,True,True,True,0.7,70-80\n',
  )
  await expect(exportCrossType(input, output, cells)).rejects.toThrow('duplicate')
})

it('accepts legacy non-finite values in unrelated analysis metadata', async () => {
  const { input, output, combined, cells } = await fixture()
  const path = join(combined, 'analysis_metadata.json')
  const source = await readFile(path, 'utf8')
  await writeFile(path, source.slice(0, -1) + ',"unused_residual_metric":NaN}')
  await expect(exportCrossType(input, output, cells)).resolves.toBeUndefined()
})
