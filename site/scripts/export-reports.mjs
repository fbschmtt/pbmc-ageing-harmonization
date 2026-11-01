import { createHash } from 'node:crypto'
import { execFileSync } from 'node:child_process'
import { fileURLToPath } from 'node:url'
import { access, mkdir, readFile, readdir, writeFile } from 'node:fs/promises'
import { dirname, join, relative, resolve, isAbsolute } from 'node:path'

async function optional(path) {
  try {
    return await readFile(path)
  } catch (error) {
    if (error.code === 'ENOENT') return null
    throw error
  }
}

export async function studyReportNames(input) {
  const directory = resolve(input, '../../harmonization/harmonization_qc')
  let entries
  try {
    entries = await readdir(directory, { withFileTypes: true })
  } catch (error) {
    if (error.code === 'ENOENT') return []
    throw error
  }
  return entries
    .filter((entry) => entry.isDirectory() && entry.name !== 'merged')
    .map((entry) => `../../harmonization/harmonization_qc/${entry.name}/report.html`)
    .sort()
}

// HTML is stored under a release-specific path so cached reports cannot be
// mistaken for another release. Companion downloads are intentionally omitted.
export async function exportReports(input, stage, releaseId, cellTypes) {
  const reports = []
  const candidates = [
    ...(await studyReportNames(input)).map((name) => ({
      source: resolve(input, name),
      name: name.replace('../../', ''),
      label: `${name.split('/').at(-2)} harmonization diagnostics`,
      study: name.split('/').at(-2),
    })),
    ...cellTypes.map((cell) => ({
      source: join(input, cell.id, 'differential_expression', 'report.html'),
      name: `analysis/per_type/${cell.id}/differential_expression/report.html`,
      label: `${cell.label} differential expression`,
      cellType: cell.label,
      warning:
        cell.warning ||
        (cell.analysisMode === 'test_only' ? 'Test-only results are not biological evidence.' : null),
    })),
    {
      source: join(input, '..', 'combined', 'trajectory', 'report.html'),
      name: 'analysis/combined/trajectory/report.html',
      label: 'Cross-cell-type trajectories',
      cellType: null,
      warning:
        cellTypes
          .map(
            (cell) =>
              cell.warning ||
              (cell.analysisMode === 'test_only' ? 'Test-only results are not biological evidence.' : null),
          )
          .find(Boolean) || null,
    },
  ]
  for (const candidate of candidates) {
    const bytes = await optional(candidate.source)
    if (!bytes) continue
    const path = `reports/${releaseId}/${candidate.name}`
    const target = join(stage, path)
    await mkdir(dirname(target), { recursive: true })
    await writeFile(target, bytes)
    reports.push({
      path,
      kind: 'report',
      label: candidate.label,
      cellType: candidate.cellType,
      study: candidate.study,
      warning: candidate.warning,
      bytes: bytes.length,
      sha256: createHash('sha256').update(bytes).digest('hex'),
    })
  }
  for (const report of reports) {
    const html = (await readFile(join(stage, report.path))).toString('utf8')
    const missing = new Set()
    for (const match of html.matchAll(/<a\b[^>]*\bhref\s*=\s*(["'])(.*?)\1/gi)) {
      const href = match[2]
      if (!href || /^(?:[a-z][a-z0-9+.-]*:|#|\/\/)/i.test(href)) continue
      const name = href.split(/[?#]/, 1)[0]
      if (!name) continue
      const target = resolve(stage, dirname(report.path), name)
      const rel = relative(stage, target)
      if (rel.startsWith('..') || isAbsolute(rel)) {
        missing.add(href)
        continue
      }
      try {
        await access(target)
      } catch (error) {
        if (error.code !== 'ENOENT') throw error
        missing.add(href)
      }
    }
    report.unavailableCompanions = [...missing].sort()
  }
  if (reports.length) {
    const reportRoot = join(stage, 'reports', releaseId)
    const script = fileURLToPath(new URL('../../scripts/render_report_index.py', import.meta.url))
    execFileSync('python3', [script, '--outdir', reportRoot], { stdio: 'inherit' })
    const indexPath = `reports/${releaseId}/index.html`
    const bytes = await readFile(join(stage, indexPath))
    reports.push({
      path: indexPath,
      kind: 'report_index',
      label: 'Analysis report index',
      bytes: bytes.length,
      sha256: createHash('sha256').update(bytes).digest('hex'),
    })
  }
  return reports
}
