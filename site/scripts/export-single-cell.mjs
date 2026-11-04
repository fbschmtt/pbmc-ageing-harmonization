import { mkdir, readFile, readdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'

async function optional(path) {
  try { return await readFile(path, 'utf8') } catch (error) {
    if (error.code === 'ENOENT') return null
    throw error
  }
}

export async function exportSingleCell(input, stage) {
  const entries = await readdir(input, { withFileTypes: true })
  const manifestText = await optional(join(input, 'cell_type_manifest.json'))
  let selected = null
  if (manifestText) {
    const manifest = JSON.parse(manifestText)
    if (!Array.isArray(manifest.cell_types) || !manifest.cell_types.every((item) =>
      typeof item.slug === 'string' && typeof item.cell_type === 'string' && typeof item.status === 'string'))
      throw new Error('Invalid single-cell analysis manifest')
    selected = new Map(manifest.cell_types.filter((item) => item.status === 'complete')
      .map((item) => [item.slug, item.cell_type]))
  }
  const cellTypes = []
  for (const entry of entries) {
    if (!entry.isDirectory() || !/^[a-z0-9]+(?:-[a-z0-9]+)*$/.test(entry.name)) continue
    if (selected && !selected.has(entry.name)) continue
    const directory = join(input, entry.name, 'analysis')
    const [reportText, notebookText, markerText] = await Promise.all([
      optional(join(directory, 'report.json')),
      optional(join(directory, 'executed.ipynb')),
      optional(join(directory, 'markers.tsv')),
    ])
    if (!reportText || !notebookText || !markerText) continue
    const report = JSON.parse(reportText)
    if (report.status !== 'complete' || typeof report.cell_type !== 'string') continue
    if (selected && selected.get(entry.name) !== report.cell_type)
      throw new Error(`Single-cell manifest label mismatch: ${entry.name}`)
    const notebook = JSON.parse(notebookText)
    const cells = notebook.cells || []
    const heading = cells.findIndex((cell) => cell.cell_type === 'markdown' &&
      (Array.isArray(cell.source) ? cell.source.join('') : cell.source || '').includes('### Harmony-Integrated UMAP'))
    const figure = heading < 0 ? null : cells[heading + 1]?.outputs?.find((output) => output.data?.['image/png'])
    const image = figure?.data?.['image/png']
    if (!image) continue
    const encoded = Array.isArray(image) ? image.join('') : image
    if (!Buffer.from(encoded, 'base64').subarray(0, 8).equals(Buffer.from('89504e470d0a1a0a', 'hex'))) continue
    const fractionCell = cells.find((cell) => (Array.isArray(cell.source) ? cell.source.join('') : cell.source || '')
      .includes('cluster_fractions = sample_cluster_fractions('))
    const fractionHasStudyTrends = cells.some((cell) => (Array.isArray(cell.source) ? cell.source.join('') : cell.source || '')
      .includes('# study-specific linear trends'))
    const fractionIncludesZeros = cells.some((cell) => (Array.isArray(cell.source) ? cell.source.join('') : cell.source || '')
      .includes('# cluster-fraction zeros included'))
    let fractionImage = null
    let fractionNote = null
    let fractionDetails = null
    for (const output of fractionCell?.outputs || []) {
      const details = output.data?.['application/vnd.pbmc.cluster-fractions+json']
      if (details && Array.isArray(details.clusters)) fractionDetails = details.clusters
      const html = output.data?.['text/html']
      const text = Array.isArray(html) ? html.join('') : html || ''
      if (text.includes('Show cluster-composition curves')) {
        const match = text.match(/src=["']data:image\/png;base64,([A-Za-z0-9+/=\s]+)["']/)
        if (match) fractionImage = match[1].replace(/\s/g, '')
      }
      const markdown = output.data?.['text/markdown']
      const message = Array.isArray(markdown) ? markdown.join('') : markdown || ''
      if (message.includes('excluded from cluster fractions')) fractionNote = message.replace(/^_|_$/g, '')
    }
    // New reports provide one image per cluster; keep the combined figure only for older reports.
    if (fractionDetails !== null) fractionImage = null
    const lines = markerText.trim().split(/\r?\n/)
    const headers = lines.shift().split('\t')
    const column = (name) => headers.indexOf(name)
    if (['names', 'logfoldchanges', 'pvals_adj', 'group'].some((name) => column(name) < 0))
      throw new Error(`Invalid marker table for ${entry.name}`)
    const clusters = new Map()
    for (const line of lines) {
      const fields = line.split('\t')
      const gene = fields[column('names')], cluster = fields[column('group')], raw = fields[column('logfoldchanges')]
      const effect = raw ? Number(raw) : NaN
      const pText = fields[column('pvals_adj')], padj = pText ? Number(pText) : null
      if (!gene || !cluster) continue
      if (!clusters.has(cluster)) clusters.set(cluster, [])
      if (!Number.isFinite(effect) || effect === 0 ||
          (padj !== null && (!Number.isFinite(padj) || padj < 0 || padj > 1))) continue
      clusters.get(cluster).push({ gene, effect, padj })
    }
    const groups = [...clusters].sort(([a], [b]) => a.localeCompare(b, undefined, { numeric: true }))
      .map(([cluster, markers]) => ({ cluster,
        higher: markers.filter((item) => item.effect > 0).sort((a, b) => b.effect - a.effect).slice(0, 8),
        lower: markers.filter((item) => item.effect < 0).sort((a, b) => a.effect - b.effect).slice(0, 8),
      }))
    const path = `single-cell/${entry.name}.json`
    await mkdir(join(stage, 'single-cell'), { recursive: true })
    await writeFile(join(stage, path), JSON.stringify({ schemaVersion: 1, cellType: report.cell_type,
      nCells: report.n_cells, nClusters: report.n_clusters, integration: report.embedding?.integration?.method || 'none',
      image: encoded, fractionImage, fractionNote, fractionHasStudyTrends, fractionIncludesZeros, fractionDetails, groups }))
    cellTypes.push({ id: entry.name, label: report.cell_type, path })
  }
  if (cellTypes.length) await writeFile(join(stage, 'single-cell-index.json'), JSON.stringify({ schemaVersion: 1, cellTypes }))
}
