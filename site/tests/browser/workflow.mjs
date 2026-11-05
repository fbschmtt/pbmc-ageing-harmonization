import { test } from 'node:test'
import assert from 'node:assert/strict'
import { createServer } from 'node:http'
import { mkdtemp, mkdir, readFile, rm } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { dirname, extname, join, resolve, sep } from 'node:path'
import { fileURLToPath } from 'node:url'
import { build } from 'vite'
import { chromium } from '@playwright/test'
import { createRelease } from '../../scripts/release.mjs'
import { installOfflineRoutes } from './network.mjs'

const site = resolve(dirname(fileURLToPath(import.meta.url)), '../..')
const input = resolve(process.env.SITE_TEST_INPUT || join(site, '../output/test/analysis/per_type'))
const runRoot = resolve(input, '../..')

async function choosePageFromGroup(page, parent, child) {
  const menu = page.locator('.site-nav .nav-menu').filter({ hasText: parent })
  await menu.locator('summary').click()
  // The submenu button's accessible name includes its descriptive second line.
  await menu.locator('.nav-menu-panel button').filter({ hasText: child }).click()
}

test(
  'workflow outputs produce a browsable atlas with harmonization and DE',
  { timeout: 180000 },
  async () => {
    const run = JSON.parse(await readFile(join(runRoot, 'run_manifest.json'), 'utf8'))
    assert.ok(
      Array.isArray(run.selected_studies) && run.selected_studies.length,
      'No studies selected in test run manifest',
    )
    const root = await mkdtemp(join(tmpdir(), 'atlas-workflow-browser-'))
    let server, browser
    try {
      const publicDir = join(root, 'public'),
        dist = join(root, 'dist')
      await mkdir(publicDir, { recursive: true })
      const release = await createRelease(input, join(publicDir, 'data'))
      assert.ok(
        release.inventory.some((cell) => cell.syntheticTestData === true),
        'Expected synthetic DE models',
      )
      for (const study of run.selected_studies) {
        assert.ok(
          release.context.studies.some((item) => item.id === study && item.harmonizationAvailable),
          `Missing study context: ${study}`,
        )
        assert.ok(
          release.assets.some((item) => item.study === study && item.kind === 'report'),
          `Missing study report: ${study}`,
        )
      }
      const figures = release.assets.filter(
        (item) =>
          ['study_umap', 'study_heatmap'].includes(item.kind) && run.selected_studies.includes(item.study),
      )
      assert.ok(figures.length > 0, 'No study diagnostic figures exported from workflow notebooks')
      const localIndex = JSON.parse(await readFile(join(publicDir, 'data/single-cell-index.json'), 'utf8'))
      assert.ok(localIndex.cellTypes.some((item) => item.label === 'CD14 monocyte'), 'Missing single-cell summary')
      process.env.ATLAS_TEST_PUBLIC_DIR = publicDir
      await build({
        configFile: join(site, 'tests/browser/vite.config.mjs'),
        build: { outDir: dist },
        logLevel: 'error',
      })
      const mime = {
        '.html': 'text/html',
        '.js': 'text/javascript',
        '.css': 'text/css',
        '.json': 'application/json',
        '.png': 'image/png',
        '.svg': 'image/svg+xml',
      }
      server = createServer(async (req, res) => {
        try {
          const name = decodeURIComponent(new URL(req.url, 'http://localhost').pathname)
          const path = resolve(dist, `.${name === '/' ? '/index.html' : name}`)
          if (!path.startsWith(dist + sep)) throw new Error('Invalid path')
          const bytes = await readFile(path)
          res.writeHead(200, { 'content-type': mime[extname(path)] || 'application/octet-stream' })
          res.end(bytes)
        } catch {
          res.writeHead(404)
          res.end('not found')
        }
      })
      await new Promise((resolve, reject) => {
        server.once('error', reject)
        server.listen(0, '127.0.0.1', resolve)
      })
      browser = await chromium.launch({
        headless: true,
        executablePath: process.env.SITE_BROWSER_EXECUTABLE || undefined,
      })
      const page = await browser.newPage()
      const appUrl = `http://127.0.0.1:${server.address().port}/`
      await installOfflineRoutes(page, appUrl)
      await page.goto(appUrl)
      await page.getByRole('button', { name: 'Datasets', exact: true }).click()
      for (const study of run.selected_studies) {
        const index = release.context.studies.findIndex((item) => item.id === study)
        assert.notEqual(index, -1, `Missing dataset card: ${study}`)
        const link = page.locator(`.study-index a[href="#study-section-${index}"]`)
        await link.waitFor({ state: 'visible', timeout: 10000 })
        assert.equal((await link.textContent()).trim(), study)
        const card = page.locator(`#study-section-${index}`)
        await card.waitFor({ state: 'visible', timeout: 10000 })
        await card.locator('h3').waitFor({ state: 'visible', timeout: 10000 })
        assert.equal(
          (await card.locator('h3').textContent()).trim(),
          release.context.studies[index].publicationTitle || study,
        )
      }
      for (const figure of figures) {
        const displayed = page.locator(`img[src$="${figure.path}"]`)
        await displayed.scrollIntoViewIfNeeded()
        await displayed.waitFor({ state: 'visible' })
        await page.waitForFunction(
          (src) =>
            [...document.images].some((img) => img.src.endsWith(src) && img.complete && img.naturalWidth > 0),
          figure.path,
        )
      }
      await choosePageFromGroup(page, 'Differential expression', 'By cell type')
      await page.locator('.de-warning').waitFor({ state: 'visible', timeout: 10000 })
      await page.locator('.de-statistics').waitFor({ state: 'visible', timeout: 10000 })
      await page.getByRole('button', { name: 'Cell states', exact: true }).click()
      await page.getByRole('img', { name: /^UMAP of / }).waitFor({ state: 'visible', timeout: 10000 })
      assert.ok(await page.locator('.single-cell-stats strong').count() >= 2)
    } finally {
      if (browser) await browser.close()
      if (server?.listening)
        await new Promise((resolve, reject) => server.close((error) => (error ? reject(error) : resolve())))
      await rm(root, { recursive: true, force: true })
    }
  },
)
