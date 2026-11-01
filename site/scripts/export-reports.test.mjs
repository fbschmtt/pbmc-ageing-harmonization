import { afterEach, describe, expect, it } from 'vitest'
import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { exportReports } from './export-reports.mjs'

let root
afterEach(async () => {
  if (root) await rm(root, { recursive: true, force: true })
})
describe('HTML diagnostics', () => {
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
