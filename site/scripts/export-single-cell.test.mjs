import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { afterEach, expect, it } from 'vitest'
import { exportSingleCell } from './export-single-cell.mjs'
import { writeSingleCellAssets, png } from '../tests/single-cell-fixture.mjs'

const roots = []
afterEach(async () => { await Promise.all(roots.splice(0).map((root) => rm(root, { recursive: true, force: true }))) })

async function analysis(input, id, label, assets = true) {
  const path = join(input, id, 'analysis')
  await mkdir(path, { recursive: true })
  await writeFile(join(path, 'report.json'), JSON.stringify({ status: 'complete', cell_type: label, n_cells: 20, n_clusters: 1 }))
  await writeFile(join(path, 'markers.tsv'), 'group\tnames\tlogfoldchanges\tpvals_adj\n0\tUP\t2\t0.01\n0\tDOWN\t-3\t0.02\n')
  if (assets) await writeSingleCellAssets(path, label)
}

it('exports explicit assets without needing notebook headings, cells or source comments', async () => {
  const root = await mkdtemp(join(tmpdir(), 'single-cell-export-')); roots.push(root)
  const input = join(root, 'per_type'), stage = join(root, 'stage')
  await mkdir(stage)
  await analysis(input, 'selected', 'Selected')
  await analysis(input, 'stale', 'Stale')
  await writeFile(join(input, 'cell_type_manifest.json'), JSON.stringify({ cell_types: [
    { slug: 'selected', cell_type: 'Selected', status: 'complete' },
  ] }))
  await exportSingleCell(input, stage)
  const index = JSON.parse(await readFile(join(stage, 'single-cell-index.json'), 'utf8'))
  expect(index.cellTypes.map((item) => item.id)).toEqual(['selected'])
  const payload = JSON.parse(await readFile(join(stage, 'single-cell/selected.json'), 'utf8'))
  expect(payload.groups[0].higher[0]).toMatchObject({ gene: 'UP', effect: 2 })
  expect(payload.groups[0].lower[0]).toMatchObject({ gene: 'DOWN', effect: -3 })
  expect(payload.fractionDetails).toBeNull()
  expect(payload.image).toMatch(/^single-cell-images\/[a-f0-9]{64}\.png$/)
  expect(await readFile(join(stage, payload.image))).toEqual(Buffer.from(png, 'base64'))
})

it('requires rerendering older reports without explicit site assets', async () => {
  const root = await mkdtemp(join(tmpdir(), 'single-cell-export-')); roots.push(root)
  const input = join(root, 'per_type'), stage = join(root, 'stage')
  await mkdir(stage)
  await analysis(input, 'selected', 'Selected', false)
  await expect(exportSingleCell(input, stage)).rejects.toThrow('rerender this cell-type report')
})

it('rejects changed images even when report metadata is unchanged', async () => {
  const root = await mkdtemp(join(tmpdir(), 'single-cell-export-')); roots.push(root)
  const input = join(root, 'per_type'), stage = join(root, 'stage')
  await mkdir(stage)
  await analysis(input, 'selected', 'Selected')
  const directory = join(input, 'selected/analysis')
  const manifest = JSON.parse(await readFile(join(directory, 'site-assets.json'), 'utf8'))
  await writeFile(join(directory, manifest.image), 'corrupt image')
  await expect(exportSingleCell(input, stage)).rejects.toThrow('checksum mismatch')
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
