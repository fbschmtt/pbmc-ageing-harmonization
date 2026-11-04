import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { afterEach, expect, it } from 'vitest'
import { exportSingleCell } from './export-single-cell.mjs'

const roots = []
const png = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+/a1sAAAAASUVORK5CYII='
afterEach(async () => { await Promise.all(roots.splice(0).map((root) => rm(root, { recursive: true, force: true }))) })

async function analysis(input, id, label, notebook = true) {
  const path = join(input, id, 'analysis')
  await mkdir(path, { recursive: true })
  await writeFile(join(path, 'report.json'), JSON.stringify({ status: 'complete', cell_type: label, n_cells: 20, n_clusters: 1 }))
  await writeFile(join(path, 'markers.tsv'), 'group\tnames\tlogfoldchanges\tpvals_adj\n0\tUP\t2\t0.01\n0\tDOWN\t-3\t0.02\n')
  if (notebook) await writeFile(join(path, 'executed.ipynb'), JSON.stringify({ cells: [
    { cell_type: 'markdown', source: ['### Harmony-Integrated UMAP'] },
    { cell_type: 'code', outputs: [{ data: { 'image/png': png } }] },
  ] }))
}

it('exports only complete manifest entries and preserves older reports without fraction plots', async () => {
  const root = await mkdtemp(join(tmpdir(), 'single-cell-export-')); roots.push(root)
  const input = join(root, 'per_type'), stage = join(root, 'stage')
  await mkdir(stage)
  await analysis(input, 'selected', 'Selected')
  await analysis(input, 'stale', 'Stale')
  await analysis(input, 'missing', 'Missing', false)
  await writeFile(join(input, 'cell_type_manifest.json'), JSON.stringify({ cell_types: [
    { slug: 'selected', cell_type: 'Selected', status: 'complete' },
    { slug: 'missing', cell_type: 'Missing', status: 'complete' },
  ] }))
  await exportSingleCell(input, stage)
  const index = JSON.parse(await readFile(join(stage, 'single-cell-index.json'), 'utf8'))
  expect(index.cellTypes.map((item) => item.id)).toEqual(['selected'])
  const payload = JSON.parse(await readFile(join(stage, 'single-cell/selected.json'), 'utf8'))
  expect(payload.groups[0].higher[0]).toMatchObject({ gene: 'UP', effect: 2 })
  expect(payload.groups[0].lower[0]).toMatchObject({ gene: 'DOWN', effect: -3 })
  expect(payload.fractionDetails).toBeNull()
})

it('rejects a report that disagrees with its manifest', async () => {
  const root = await mkdtemp(join(tmpdir(), 'single-cell-export-')); roots.push(root)
  const input = join(root, 'per_type'), stage = join(root, 'stage')
  await mkdir(stage)
  await analysis(input, 'selected', 'Actual')
  await writeFile(join(input, 'cell_type_manifest.json'), JSON.stringify({ cell_types: [
    { slug: 'selected', cell_type: 'Different', status: 'complete' },
  ] }))
  await expect(exportSingleCell(input, stage)).rejects.toThrow('manifest label mismatch')
})
