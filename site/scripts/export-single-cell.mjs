import { mkdir, readFile, readdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'
import { readSingleCellAssets } from './single-cell-assets.mjs'
import { writeSingleCellImage } from './png.mjs'

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
    const [reportText, markerText] = await Promise.all([
      optional(join(directory, 'report.json')),
      optional(join(directory, 'markers.tsv')),
    ])
    if (!reportText || !markerText) continue
    const report = JSON.parse(reportText)
    if (report.status !== 'complete' || typeof report.cell_type !== 'string') continue
    if (selected && selected.get(entry.name) !== report.cell_type)
      throw new Error(`Single-cell manifest label mismatch: ${entry.name}`)
    const assets = await readSingleCellAssets(directory, report.cell_type)
    if (!assets.image) continue
    const encoded = await readFile(join(directory, assets.image))
    let ageDeImage = assets.ageDeImage ? await readFile(join(directory, assets.ageDeImage)) : null
    const { ageDeNote, fractionNote, fractionHasStudyTrends, fractionIncludesZeros } = assets
    let fractionDetails = assets.fractionDetails
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
        higher: markers.filter((item) => item.effect > 0).sort((a, b) => b.effect - a.effect).slice(0, 10),
        lower: markers.filter((item) => item.effect < 0).sort((a, b) => a.effect - b.effect).slice(0, 10),
      }))
    const imageAsset = await writeSingleCellImage(stage, encoded.toString('base64'))
    if (ageDeImage) ageDeImage = await writeSingleCellImage(stage, ageDeImage.toString('base64'))
    if (fractionDetails) fractionDetails = await Promise.all(fractionDetails.map(async (detail) => ({
      ...detail, image: await writeSingleCellImage(stage, (await readFile(join(directory, detail.image))).toString('base64')),
    })))
    const path = `single-cell/${entry.name}.json`
    await mkdir(join(stage, 'single-cell'), { recursive: true })
    await writeFile(join(stage, path), JSON.stringify({ schemaVersion: 1, cellType: report.cell_type,
      nCells: report.n_cells, nClusters: report.n_clusters, integration: report.embedding?.integration?.method || 'none',
      image: imageAsset, ageDeImage, ageDeNote, fractionNote, fractionHasStudyTrends, fractionIncludesZeros, fractionDetails, groups }))
    cellTypes.push({ id: entry.name, label: report.cell_type, path })
  }
  if (cellTypes.length) await writeFile(join(stage, 'single-cell-index.json'), JSON.stringify({ schemaVersion: 1, cellTypes }))
}
