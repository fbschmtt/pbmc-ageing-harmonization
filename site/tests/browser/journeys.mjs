import { test } from 'node:test'
import assert from 'node:assert/strict'
import { build } from 'vite'
import { createServer } from 'node:http'
import { mkdtemp, readFile, rm } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { dirname, extname, join, resolve, sep } from 'node:path'
import { fileURLToPath } from 'node:url'
import { gzipSync } from 'node:zlib'
import { chromium } from '@playwright/test'
import { makeFixture } from './fixture.mjs'
import { verifyRelease } from '../../scripts/verify-release.mjs'

const site = resolve(dirname(fileURLToPath(import.meta.url)), '../..')
const types = {
  '.html': 'text/html',
  '.js': 'text/javascript',
  '.css': 'text/css',
  '.json': 'application/json',
  '.svg': 'image/svg+xml',
}

async function serve(root) {
  const server = createServer(async (req, res) => {
    try {
      const pathname = decodeURIComponent(new URL(req.url, 'http://localhost').pathname)
      const path = resolve(root, `.${pathname === '/' ? '/index.html' : pathname}`)
      if (!path.startsWith(root + sep)) throw new Error('invalid path')
      const data = await readFile(path)
      res.writeHead(200, { 'content-type': types[extname(path)] || 'application/octet-stream' })
      res.end(data)
    } catch {
      res.writeHead(404)
      res.end('not found')
    }
  })
  await new Promise((resolve, reject) => {
    server.once('error', reject)
    server.listen(0, '127.0.0.1', resolve)
  })
  return { server, url: `http://127.0.0.1:${server.address().port}/` }
}
const visible = async (locator) => locator.waitFor({ state: 'visible', timeout: 10000 })
const menuTrigger = (page, label) =>
  page
    .locator('.choice-field')
    .filter({ has: page.locator('.choice-label', { hasText: label }) })
    .locator('summary')
const chooseMenu = async (page, label, option) => {
  await menuTrigger(page, label).click()
  await page
    .getByRole('menu', { name: label, exact: true })
    .getByRole('menuitemradio', { name: option, exact: true })
    .click()
}
const assertMenu = async (page, label, value) => {
  assert.equal(await menuTrigger(page, label).getAttribute('aria-label'), `${label}: ${value}`)
}

const group = process.env.SITE_BROWSER_GROUP || 'all'
const filter = (process.env.SITE_BROWSER_FILTER || '').toLowerCase()
if (!['all', 'navigation', 'release'].includes(group)) throw new Error(`Unknown browser test group: ${group}`)

// One real exported and built release is shared across independent browser contexts.
test('built synthetic release supports atlas journeys and fails safely', { timeout: 120000 }, async (t) => {
  const root = await mkdtemp(join(tmpdir(), 'atlas-browser-'))
  let server, browser
  let matched = 0
  const runJourney = async (testGroup, name, fn) => {
    if ((group !== 'all' && group !== testGroup) || !name.toLowerCase().includes(filter)) return
    matched++
    await t.test(name, fn)
  }
  try {
    const { publicDir, manifest } = await makeFixture(root)
    const dist = join(root, 'dist')
    process.env.ATLAS_TEST_PUBLIC_DIR = publicDir
    await build({
      configFile: join(site, 'tests/browser/vite.config.mjs'),
      build: { outDir: dist },
      logLevel: 'error',
    })
    const built = await verifyRelease(join(dist, 'data'))
    assert.equal(built.releaseId, manifest.releaseId)
    const report = built.assets.find((asset) => asset.kind === 'report' && asset.cellType === 'CD14 monocyte')
    assert.ok(report)
    assert.equal(
      built.assets.some((asset) => asset.kind === 'metadata'),
      false,
    )
    assert.match(await readFile(join(dist, 'data', report.path), 'utf8'), /Fixture DE report/)
    const served = await serve(dist)
    server = served.server
    browser = await chromium.launch({
      headless: true,
      executablePath: process.env.SITE_BROWSER_EXECUTABLE || undefined,
    })

    await runJourney('release', 'fetches only the selected DE model and uses compressed assets', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        const requests = []
        page.on('request', (request) => requests.push(new URL(request.url()).pathname))
        await page.goto(served.url)
        assert.equal(
          requests.some((path) => path.includes('/data/de-models/')),
          false,
        )
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(page.locator('.de-statistics'))
        const catalogue = JSON.parse(await readFile(join(dist, 'data/de-index.json'), 'utf8'))
        const cell = catalogue.cellTypes.find((item) => item.label === 'CD14 monocyte')
        const other = cell.models.find((model) => model.model === 'per_study')
        const otherPaths = [other.recordsAsset, ...Object.values(other.recordBuckets || {})]
        await page.getByLabel('Highlight genes').fill('GENE3')
        await page.waitForTimeout(1400)
        const wasRequested = (asset) =>
          requests.some((path) => path.endsWith(`/${asset}`) || path.endsWith(`/${asset}.gz.bin`))
        assert.equal(otherPaths.some(wasRequested), false, 'Unselected model was prefetched')
        assert.ok(requests.some((path) => path.includes('/data/de-models/') && path.endsWith('.gz.bin')))
        await chooseMenu(page, 'Model and contrast', 'study_a · Age (per decade)')
        await visible(page.getByText('Effect: age · per decade.', { exact: false }))
        assert.equal(wasRequested(other.recordsAsset), true, 'Selected model was not fetched')
        await visible(page.getByText('1 of 1 studies significant', { exact: true }))
        assert.ok((await page.locator('svg[aria-label^="Volcano plot"] title').allTextContents())
          .some((title) => title.includes('1 of 1 studies significant')))
      } finally {
        await context.close()
      }
    })

    await runJourney('release', 'filters legacy inline results and explains their model-loading limitation', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        const legacy = JSON.parse(await readFile(join(dist, 'data/de-index.json'), 'utf8'))
        for (const cell of legacy.cellTypes)
          for (const model of cell.models) {
            const paths = [model.recordsAsset, ...Object.values(model.recordBuckets || {})]
            model.records = []
            for (const path of paths) {
              const chunk = JSON.parse(await readFile(join(dist, 'data', path), 'utf8'))
              model.records.push(...chunk.cellTypes[0].models[0].records)
            }
            delete model.recordsAsset
            delete model.recordBuckets
          }
        await page.route('**/data/de-index.json*', (route) =>
          route.request().url().endsWith('.gz.bin')
            ? route.fulfill({
                body: gzipSync(JSON.stringify(legacy)),
                contentType: 'application/octet-stream',
              })
            : route.fulfill({ json: legacy }),
        )
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(
          page.getByText('This older export contains all models in one download.', { exact: false }),
        )
        const plot = page.locator('svg[aria-label^="Volcano plot"]')
        await visible(plot)
        assert.equal(await plot.locator('circle').count(), 2)
        assert.equal(await plot.locator('title').filter({ hasText: 'GENE3' }).count(), 0)
        await page.getByLabel('Find a gene').fill('GENE3')
        await chooseMenu(page, 'Matching genes (1)', 'GENE3')
        assert.equal(await plot.locator('title').filter({ hasText: 'GENE3' }).count(), 1)
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'retrieves a non-significant gene on explicit search', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(page.locator('.de-statistics'))
        await page.getByLabel('Find a gene').fill('GENE3')
        await chooseMenu(page, 'Matching genes (1)', 'GENE3')
        await visible(page.getByRole('heading', { name: 'GENE3', exact: true }))
        await visible(page.locator('.de-statistics dd').first().filter({ hasText: '0.230' }))
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'retrieves non-significant retained trajectories for an explicit gene', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByText('Age-bin trajectories', { exact: false }).first().click()
        await page.getByRole('button', { name: /Genes.*Compare a gene/ }).click()
        await visible(page.getByLabel('Gene symbol'))
        await page.getByLabel('Gene symbol').fill('GENE2')
        await page
          .locator('.profile-scatter title')
          .filter({ hasText: 'Naive CD4 T cell' })
          .waitFor({ state: 'attached' })
        assert.equal(await page.locator('.profile-scatter circle').filter({ has: page.locator('title', { hasText: 'Naive CD4 T cell' }) }).getAttribute('fill'), '#597FC6')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'overview explains the workflow on desktop and mobile', async () => {
      const context = await browser.newContext({ viewport: { width: 375, height: 800 } })
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        const workflow = page.getByRole('region', { name: 'How the atlas comes together' })
        await visible(workflow)
        assert.equal(await workflow.locator('li').count(), 4)
        await visible(workflow.getByText('Sum gene counts by sample and cell type; count cells for composition.'))
        assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), true)
        await workflow.getByRole('link', { name: 'Composition', exact: true }).click()
        await visible(page.locator('.pipeline-composition'))
        await page.getByRole('button', { name: 'Overview', exact: true }).click()
        await page.setViewportSize({ width: 1280, height: 800 })
        await visible(workflow)
        assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), true)
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'single-cell landscape shows embedding and ranked markers', async () => {
      const context = await browser.newContext({ viewport: { width: 390, height: 844 } })
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Single-cell', exact: true }).click()
        await visible(page.getByRole('heading', { name: 'Single-cell landscape' }))
        await visible(page.getByRole('img', { name: 'UMAP of CD14 monocyte, colored by local cluster' }))
        await visible(page.getByText('TESTMARKER', { exact: true }))
        await visible(page.getByText('OTHER', { exact: true }))
        await visible(page.getByRole('img', { name: /Cluster 0 sample-level fractions across age/ }))
        await visible(page.getByRole('cell', { name: '25.0%' }))
        await visible(page.getByRole('cell', { name: '+1.50' }))
        await visible(page.getByText('log₂ FC').first())
        await visible(page.getByText(/Dashed lines are per-study trends/))
        await visible(page.getByText(/Eligible samples with no cells in a cluster are included at zero/))
        await visible(page.getByText(/3 sample\(s\) excluded/))
        await page.evaluate(() => document.fonts.ready)
        const layout = () => page.evaluate(() => {
          const fraction = document.querySelector('.single-cell-fractions')
          const frame = document.querySelector('.single-cell-plot-frame')
          const table = document.querySelector('.single-cell-fit-table')
          return {
            fractionTop: fraction.getBoundingClientRect().top + window.scrollY,
            frameHeight: frame.getBoundingClientRect().height,
            tableHeight: table.getBoundingClientRect().height,
          }
        })
        const before = await layout()
        await chooseMenu(page, 'Cluster', 'Cluster 1')
        await visible(page.getByText('SECOND', { exact: true }))
        await visible(page.getByRole('img', { name: /Cluster 1 sample-level fractions across age/ }))
        await visible(page.getByRole('cell', { name: '75.0%' }))
        const after = await layout()
        assert.ok(Math.abs(before.fractionTop - after.fractionTop) < 1, 'Fraction panel moved')
        assert.ok(Math.abs(before.frameHeight - after.frameHeight) < 1, 'Plot frame changed height')
        assert.ok(Math.abs(before.tableHeight - after.tableHeight) < 1, 'Fit table changed height')
        assert.equal(await page.getByText('TESTMARKER', { exact: true }).count(), 0)
        assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), true)
      } finally { await context.close() }
    })

    await runJourney('navigation', 'age model guide explains the three models on mobile', async () => {
      const context = await browser.newContext({ viewport: { width: 390, height: 844 } })
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Modelling', exact: true }).click()
        const guide = page.getByRole('region', { name: 'Three ways to look at age' })
        await visible(guide)
        for (const name of ['Linear in adults', 'Linear in juveniles', 'Age-bin model'])
          await visible(guide.getByRole('heading', { name, exact: true }))
        await visible(guide.getByText(/different fits/))
        assert.equal(await page.evaluate(() => document.documentElement.scrollWidth <= window.innerWidth), true)
      } finally { await context.close() }
    })

    await runJourney('navigation', 'overview cards navigate from their headings and keyboard', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('heading', { name: 'Cell composition' }).click()
        await visible(page.locator('.pipeline-composition'))
        assert.equal(new URL(page.url()).hash, '#composition')
        await page.getByRole('button', { name: 'Overview', exact: true }).click()
        const trajectory = page.locator('.entry-grid a[href="#trajectories"]')
        await trajectory.focus()
        await page.keyboard.press('Enter')
        await visible(page.locator('.explorer'))
        assert.equal(new URL(page.url()).hash, '#trajectories')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'model selection, gene handoff, composition and diagnostic report', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage(),
          errors = []
        page.on('pageerror', (error) => errors.push(error.message))
        await page.goto(served.url)
        await visible(page.getByText('SYNTHETIC TEST DATA · EXPORT', { exact: false }))
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(page.getByText('Synthetic fixture: not biological evidence.'))
        await chooseMenu(page, 'Model and contrast', 'study_a · Age (per decade)')
        await visible(page.getByText('Effect: age · per decade.', { exact: false }))
        await page.getByLabel('Find a gene').fill('GENE2')
        await chooseMenu(page, 'Matching genes (1)', 'GENE2')
        await visible(page.getByRole('heading', { name: 'GENE2', exact: true }))
        await visible(page.locator('.de-statistics dd').first().filter({ hasText: '-0.210' }))
        await page.getByText('Age-bin trajectories', { exact: false }).first().click()
        await page.getByRole('button', { name: /Genes.*Compare a gene/ }).click()
        await visible(page.locator('.profiles-page'))
        assert.equal(await page.getByLabel('Gene symbol').inputValue(), 'GENE2')
        await visible(page.getByText('metadata unavailable', { exact: false }).first())
        await page.getByText('Age-bin trajectories', { exact: false }).first().click()
        await page.getByRole('button', { name: /Cross-type trajectory map.*Explore/ }).click()
        await visible(page.locator('.pair-card').getByRole('button', { name: 'GENE2', exact: true }))
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(page.getByText('Synthetic fixture: not biological evidence.'))
        assert.equal(await page.locator('.fraction-baseline strong').textContent(), '25.00%')
        await chooseMenu(page, 'Age plot denominator', 'Within retained PBMCs')
        await visible(page.locator('.fraction-baseline strong').filter({ hasText: '40.00%' }))
        const agePlot = page.locator('svg[aria-label^="Model-adjusted absolute fractions"]')
        const tooltip = await agePlot.locator('title').allTextContents()
        assert.ok(tooltip.some((value) => value.includes('60.00% estimated fraction')))
        assert.ok((await agePlot.locator('path').count()) > 0)
        await visible(agePlot.getByText(/^40%$/))
        await page.getByRole('button', { name: 'Modelling', exact: true }).click()
        const reportLink = page.getByRole('link', { name: 'Open the report index' })
        await visible(reportLink)
        const index = built.assets.find((asset) => asset.kind === 'report_index')
        assert.ok((await reportLink.getAttribute('href')).endsWith(index.path))
        const [popup] = await Promise.all([page.waitForEvent('popup'), reportLink.click()])
        await visible(popup.getByRole('heading', { name: 'PBMC ageing analysis reports' }))
        await popup.getByRole('link', { name: 'DE & abundance' }).click()
        await visible(popup.getByRole('heading', { name: 'Fixture DE report' }))
        assert.ok(popup.url().endsWith(report.path))
        assert.deepEqual(errors, [])
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'shows release provenance, study context and sharing limitations', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Datasets', exact: true }).click()
        await visible(page.getByText('study_a', { exact: true }))
        await visible(page.getByText('Fixture annotation', { exact: false }))
        assert.equal(
          await page.getByRole('link', { name: '10.1234/fixture' }).getAttribute('href'),
          'https://doi.org/10.1234/fixture',
        )
        await visible(page.getByText('metadata comparison', { exact: false }))
        await page.getByRole('button', { name: 'Modelling', exact: true }).click()
        await visible(page.getByText('Study licences are not recorded in this release.', { exact: false }))
        await visible(
          page.getByText('Included HTML reports may contain sample-level identifiers', { exact: false }),
        )
        await page.getByRole('button', { name: 'Download', exact: true }).click()
        await visible(page.getByText('Downloads will be available in a future release.'))
        assert.equal(await page.locator('.content-page a').count(), 0)
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'supports keyboard selection and dismissal of menus', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        const trigger = menuTrigger(page, 'Cell type')
        await trigger.focus()
        await page.keyboard.press('ArrowDown')
        await page.keyboard.press('End')
        const focusedOption = await page.evaluate(() =>
          document.activeElement.textContent.trim().replace('✓', '').trim(),
        )
        await page.keyboard.press('Enter')
        await assertMenu(page, 'Cell type', focusedOption)
        assert.equal(await trigger.locator('..').evaluate((details) => details.open), false)
        await page.keyboard.press('Enter')
        await page.keyboard.press('Escape')
        assert.equal(await trigger.locator('..').evaluate((details) => details.open), false)
        assert.equal(await trigger.evaluate((element) => document.activeElement === element), true)
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'shares the selected cell type between DE and composition', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(menuTrigger(page, 'Cell type'))
        await page.getByLabel('Find a gene').fill('GENE2')
        await chooseMenu(page, 'Matching genes (1)', 'GENE2')
        await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
        assert.equal(await menuTrigger(page, 'Cell type').locator('.choice-color').evaluate((dot) => dot.style.backgroundColor), 'rgb(89, 127, 198)')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(page.locator('.pipeline-composition'))
        await assertMenu(page, 'Cell type', 'Naive CD4 T cell')
        await chooseMenu(page, 'Cell type', 'CD14 monocyte')
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(page.getByRole('heading', { name: 'GENE2', exact: true }))
        await assertMenu(page, 'Cell type', 'CD14 monocyte')
        await chooseMenu(page, 'Cell type', 'CD16 monocyte')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(
          page.getByText('Fraction summaries are unavailable for CD16 monocyte. Showing CD14 monocyte.'),
        )
        await assertMenu(page, 'Cell type', 'CD14 monocyte')
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(menuTrigger(page, 'Cell type'))
        await assertMenu(page, 'Cell type', 'CD16 monocyte')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await page.getByRole('button', { name: 'Use CD14 monocyte' }).click()
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(menuTrigger(page, 'Cell type'))
        await assertMenu(page, 'Cell type', 'CD14 monocyte')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'shares a cell type chosen in the trajectory profile', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await page.goto(served.url)
        await page.getByText('Age-bin trajectories', { exact: false }).first().click()
        await page.getByRole('button', { name: /Cell types.*Explore genes/ }).click()
        await visible(page.locator('.profiles-page'))
        await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(page.locator('.pipeline-composition'))
        await assertMenu(page, 'Cell type', 'Naive CD4 T cell')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation',
      'keeps model and denominator choices independent from cell type and navigation',
      async () => {
        const context = await browser.newContext()
        try {
          const page = await context.newPage()
          await page.goto(served.url)
          await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
          await visible(menuTrigger(page, 'Model and contrast'))
          await chooseMenu(page, 'Model and contrast', 'study_a · Age (per decade)')
          await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
          await assertMenu(page, 'Model and contrast', 'study_a · Age (per decade)')
          await page.getByLabel('Find a gene').fill('GENE2')
          await chooseMenu(page, 'Matching genes (1)', 'GENE2')
          await visible(page.locator('.de-statistics dd').first().filter({ hasText: '-0.670' }))
          await chooseMenu(page, 'Cell type', 'CD16 monocyte')
          await visible(
            page.getByText('The selected model is unavailable for CD16 monocyte.', { exact: false }),
          )
          await assertMenu(page, 'Model and contrast', 'study_a · Age (per decade) (unavailable)')
          await chooseMenu(page, 'Cell type', 'CD14 monocyte')
          await visible(page.locator('.de-statistics dd').first().filter({ hasText: '-0.210' }))
          await page.getByRole('button', { name: 'Composition', exact: true }).click()
          await visible(menuTrigger(page, 'Age plot denominator'))
          await chooseMenu(page, 'Age plot denominator', 'Within retained PBMCs')
          await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
          await assertMenu(page, 'Age plot denominator', 'Within retained PBMCs')
          await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
          await visible(menuTrigger(page, 'Model and contrast'))
          await assertMenu(page, 'Model and contrast', 'study_a · Age (per decade)')
          await assertMenu(page, 'Cell type', 'Naive CD4 T cell')
          await page.getByRole('button', { name: 'Composition', exact: true }).click()
          await visible(menuTrigger(page, 'Age plot denominator'))
          await assertMenu(page, 'Age plot denominator', 'Within retained PBMCs')
        } finally {
          await context.close()
        }
      },
    )

    await runJourney('release', 'shows loading until the DE asset arrives', async () => {
      const context = await browser.newContext()
      let resume
      const gate = new Promise((resolve) => {
        resume = resolve
      })
      try {
        const page = await context.newPage()
        await page.route('**/data/de-index.json*', async (route) => {
          await gate
          await route.continue()
        })
        await page.goto(served.url)
        await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
        await visible(page.getByRole('heading', { name: 'Loading differential expression…' }))
        resume()
        await visible(page.getByText('Synthetic fixture: not biological evidence.'))
      } finally {
        resume()
        await context.close()
      }
    })

    await runJourney('release',
      'invalid optional annotations are reported without blocking the trajectory view',
      async () => {
        const context = await browser.newContext()
        try {
          const page = await context.newPage()
          const updated = {
            ...built,
            assets: [
              ...built.assets,
              { path: 'gene-metadata.json', kind: 'metadata', bytes: 1, sha256: 'fixture' },
            ],
          }
          await page.route('**/data/release.json*', (route) =>
            route.request().url().endsWith('.gz.bin')
              ? route.fulfill({
                  body: gzipSync(JSON.stringify(updated)),
                  contentType: 'application/octet-stream',
                })
              : route.fulfill({ json: updated }),
          )
          await page.route('**/data/gene-metadata.json', (route) =>
            route.fulfill({ json: { releaseId: built.releaseId, schemaVersion: 1, genes: 'invalid' } }),
          )
          await page.goto(served.url)
          await page.getByText('Age-bin trajectories', { exact: false }).first().click()
          await page.getByRole('button', { name: /Genes.*Compare a gene/ }).click()
          await visible(page.locator('.profiles-page'))
          await visible(
            page.getByText(
              'Gene annotations could not be loaded because the release data is invalid or mismatched.',
            ),
          )
        } finally {
          await context.close()
        }
      },
    )

    for (const scenario of ['malformed', 'different release']) {
      await runJourney('release', `${scenario} DE asset shows an error instead of hanging`, async () => {
        const context = await browser.newContext()
        try {
          const page = await context.newPage()
          const errors = []
          page.on('pageerror', (error) => errors.push(error.message))
          const original = JSON.parse(await readFile(join(dist, 'data/de-index.json'), 'utf8'))
          const payload =
            scenario === 'malformed'
              ? { ...original, cellTypes: [{ ...original.cellTypes[0], models: 'broken' }] }
              : { ...original, releaseId: 'different-release' }
          await page.route('**/data/de-index.json*', (route) =>
            route.request().url().endsWith('.gz.bin')
              ? route.fulfill({
                  body: gzipSync(JSON.stringify(payload)),
                  contentType: 'application/octet-stream',
                })
              : route.fulfill({ json: payload }),
          )
          await page.goto(served.url)
          await page.getByRole('button', { name: 'Differential expression', exact: true }).click()
          await visible(
            page.getByRole('heading', {
              name:
                scenario === 'malformed'
                  ? 'Differential expression data is invalid.'
                  : 'Differential expression data belongs to a different release.',
            }),
          )
          assert.deepEqual(errors, [])
        } finally {
          await context.close()
        }
      })
    }
    if (!matched) throw new Error(`No browser journeys matched group ${group} and filter ${filter || '(none)'}`)
  } finally {
    if (browser) await browser.close()
    if (server)
      await new Promise((resolve, reject) => server.close((error) => (error ? reject(error) : resolve())))
    await rm(root, { recursive: true, force: true })
  }
})
