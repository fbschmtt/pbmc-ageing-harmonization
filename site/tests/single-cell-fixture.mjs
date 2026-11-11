import { createHash } from 'node:crypto'
import { mkdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'

export const png = 'iVBORw0KGgoAAAANSUhEUgAAAAEAAAABCAQAAAC1HAwCAAAAC0lEQVR42mP8/x8AAwMCAO+ip1sAAAAASUVORK5CYII='

export async function writeSingleCellAssets(directory, cellType, details = null) {
  const bytes = Buffer.from(png, 'base64')
  const sha256 = createHash('sha256').update(bytes).digest('hex')
  const image = `site-assets/${sha256}.png`
  await mkdir(join(directory, 'site-assets'), { recursive: true })
  await writeFile(join(directory, image), bytes)
  await writeFile(join(directory, 'site-assets.json'), JSON.stringify({
    schemaVersion: 1, cellType, image, ageDeImage: null, ageDeNote: null,
    fractionNote: details ? '3 sample(s) excluded from cluster fractions because the split cell-type denominator has fewer than 10 cells.' : null,
    fractionHasStudyTrends: details !== null, fractionIncludesZeros: details !== null,
    fractionDetails: details?.map((detail) => ({ ...detail, image })) ?? null,
    assets: { [image]: { sha256 } },
  }))
}
