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

const site = resolve(dirname(fileURLToPath(import.meta.url)), '../..')
const input = resolve(process.env.SITE_TEST_INPUT || join(site, '../output/test/analysis/per_type'))
const runRoot = resolve(input, '../..')

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
      await page.goto(`http://127.0.0.1:${server.address().port}/`)
      await page.getByRole('button', { name: 'Datasets', exact: true }).click()
      for (const study of run.selected_studies)
        await page
          .getByRole('heading', { name: study, exact: true })
          .waitFor({ state: 'visible', timeout: 10000 })
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
      await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
      await page.locator('.de-warning').waitFor({ state: 'visible', timeout: 10000 })
      await page.locator('.de-statistics').waitFor({ state: 'visible', timeout: 10000 })
      await page.getByRole('button', { name: 'Single-cell', exact: true }).click()
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
