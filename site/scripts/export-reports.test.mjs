import { afterEach, describe, expect, it } from 'vitest'
import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { crc32 } from 'node:zlib'
import { optimizePng } from './png.mjs'
import { exportReports } from './export-reports.mjs'

const png = (label) => {
  const base = Buffer.from('iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII=', 'base64')
  const text = Buffer.from(`Description\0${label}`)
  const chunk = Buffer.alloc(text.length + 12)
  chunk.writeUInt32BE(text.length)
  chunk.write('tEXt', 4)
  text.copy(chunk, 8)
  chunk.writeUInt32BE(crc32(chunk.subarray(4, -4)), chunk.length - 4)
  return Buffer.concat([base.subarray(0, -12), chunk, base.subarray(-12)])
}

let root
afterEach(async () => {
  if (root) await rm(root, { recursive: true, force: true })
})
describe('HTML diagnostics', () => {
  it('exports the dedicated AIFI L2 figure rather than the multi-panel UMAP', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-umap-'))
    const input = join(root, 'analysis/per_type')
    const qc = join(root, 'harmonization/harmonization_qc/example')
    await mkdir(input, { recursive: true })
    await mkdir(qc, { recursive: true })
    await writeFile(join(qc, 'report.html'), '<html>QC</html>')
    const plot = (bytes) => ({
      cell_type: 'code',
      outputs: [{ data: { 'image/png': png(bytes).toString('base64') } }],
    })
    await writeFile(
      join(qc, 'executed.ipynb'),
      JSON.stringify({
        cells: [
          { cell_type: 'markdown', source: ['### Stored UMAP embeddings'] },
          plot('multi-panel'),
          { cell_type: 'markdown', source: ['### AIFI L2 UMAP'] },
          plot('dedicated-l2'),
        ],
      }),
    )
    const output = join(root, 'release')
    const reports = await exportReports(input, output, 'release-a', [])
    const image = reports.find((item) => item.kind === 'study_umap')
    expect(await readFile(join(output, image.path))).toEqual(optimizePng(png('dedicated-l2')))
    await writeFile(
      join(qc, 'executed.ipynb'),
      JSON.stringify({
        cells: [{ cell_type: 'markdown', source: ['### Stored UMAP embeddings'] }, plot('multi-panel')],
      }),
    )
    const oldReports = await exportReports(input, join(root, 'old'), 'release-b', [])
    expect(oldReports.some((item) => item.kind === 'study_umap')).toBe(false)
  })
  it('discovers study reports independently of model cell types', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-studies-'))
    const input = join(root, 'analysis/per_type')
    const qc = join(root, 'harmonization/harmonization_qc/example')
    await mkdir(input, { recursive: true })
    await mkdir(qc, { recursive: true })
    await writeFile(join(qc, 'report.html'), '<html>UMAP</html>')
    const reports = await exportReports(input, join(root, 'release'), 'release-a', [])
    expect(reports.find((item) => item.study === 'example')?.label).toBe('example harmonization diagnostics')
  })
  it('ignores stale study directories outside the run manifest selection', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-selected-studies-'))
    const input = join(root, 'analysis/per_type')
    await mkdir(input, { recursive: true })
    for (const name of ['selected', 'stale']) {
      const dir = join(root, 'harmonization/harmonization_qc', name)
      await mkdir(dir, { recursive: true })
      await writeFile(join(dir, 'report.html'), '<html>QC</html>')
    }
    await writeFile(join(root, 'run_manifest.json'), JSON.stringify({ selected_studies: ['selected'] }))
    const reports = await exportReports(input, join(root, 'release'), 'release-a', [])
    expect(reports.filter((item) => item.study).map((item) => item.study)).toEqual(['selected'])
  })
  it('copies available reports under a release-specific path and preserves the cross-report link', async () => {
    root = await mkdtemp(join(tmpdir(), 'atlas-reports-'))
    const input = join(root, 'analysis/per_type'),
      dir = join(input, 'mono/differential_expression')
    const combined = join(root, 'analysis/combined/trajectory'),
      out = join(root, 'release')
    await mkdir(dir, { recursive: true })
    await mkdir(combined, { recursive: true })
    const html =
      '<html><a href="../../../combined/trajectory/report.html">trajectory</a><a href="sample-scores.csv">samples</a></html>'
    await writeFile(join(dir, 'report.html'), html)
    await writeFile(join(combined, 'report.html'), '<html>trajectory</html>')
    const reports = await exportReports(input, out, 'release-a', [
      { id: 'mono', label: 'Monocyte', analysisMode: 'test_only', warning: 'Synthetic' },
      { id: 'missing', label: 'Missing', analysisMode: 'standard', warning: null },
    ])
    expect(reports).toHaveLength(3)
    expect(reports[0].warning).toBe('Synthetic')
    expect(reports[0].unavailableCompanions).toEqual(['sample-scores.csv'])
    expect(reports[1].unavailableCompanions).toEqual([])
    expect(await readFile(join(out, reports[0].path), 'utf8')).toBe(html)
    expect(reports[1].path).toBe('reports/release-a/analysis/combined/trajectory/report.html')
    const index = reports.find((item) => item.kind === 'report_index')
    expect(index.path).toBe('reports/release-a/index.html')
    const htmlIndex = await readFile(join(out, index.path), 'utf8')
    expect(htmlIndex).toContain('analysis/per_type/mono/differential_expression/report.html')
    expect(htmlIndex).toContain('analysis/combined/trajectory/report.html')
    expect(htmlIndex).not.toContain('harmonization/harmonization_qc/merged/report.html')
  })
})
