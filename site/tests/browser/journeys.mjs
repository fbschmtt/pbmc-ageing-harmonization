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
import { installOfflineRoutes } from './network.mjs'

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
async function choosePageFromGroup(page, parent, child) {
  const menu = page.locator('.site-nav .nav-menu').filter({ hasText: parent })
  await menu.locator('summary').click()
  // Submenu buttons include a descriptive second line in their accessible name,
  // so match the visible page label within the selected menu instead.
  await menu.locator('.nav-menu-panel button').filter({ hasText: child }).click()
}
async function openApp(page, url) {
  await installOfflineRoutes(page, url)
  await page.goto(url)
}
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
const assertLegendEndpointAlignment = async (page) => {
  const bar = await page.locator('.de-support-gradient').boundingBox()
  const low = await page.locator('.de-support-endpoints span').first().boundingBox()
  const high = await page.locator('.de-support-endpoints span').last().boundingBox()
  assert.ok(bar && low && high)
  assert.ok(Math.abs(low.x - bar.x) < 2, 'Lower legend value should align to the left end of the color bar')
  assert.ok(Math.abs(high.x + high.width - bar.x - bar.width) < 2, 'Upper legend value should align to the right end of the color bar')
}
const clickPoint = async (plot, titleText) => {
  const points = plot.locator('circle')
  const index = await points.evaluateAll((circles, fragment) =>
    circles.findIndex((circle) => circle.querySelector('title')?.textContent?.includes(fragment)), titleText)
  assert.notEqual(index, -1, `No plot point title contains ${titleText}`)
  await points.nth(index).click()
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
    assert.ok(built.assets.some((asset) => asset.kind === 'trajectoryMap'))
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

    await runJourney('release', 'preloads combined age results and keeps per-study fits lazy', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        const requests = []
        page.on('request', (request) => requests.push(new URL(request.url()).pathname))
        await openApp(page, served.url)
        assert.equal(
          requests.some((path) => path.includes('/data/de-models/')),
          false,
        )
        await page.waitForTimeout(1600)
        const catalogue = JSON.parse(await readFile(join(dist, 'data/de-index.json'), 'utf8'))
        const cell = catalogue.cellTypes.find((item) => item.label === 'CD14 monocyte')
        const combined = cell.models.find((model) => model.model === 'merged' && model.covariate === 'age')
        const other = cell.models.find((model) => model.model === 'per_study')
        const otherPaths = [other.recordsAsset, ...Object.values(other.recordBuckets || {})]
        const wasRequested = (asset) =>
          requests.some((path) => path.endsWith(`/${asset}`) || path.endsWith(`/${asset}.gz.bin`))
        assert.equal(wasRequested(combined.recordsAsset), true, 'Combined age results were not preloaded')
        assert.equal(otherPaths.some(wasRequested), false, 'Per-study model was prefetched')
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(page.locator('.de-statistics'))
        const volcano = page.locator('.de-volcano-plot')
        const volcanoBox = await volcano.boundingBox()
        assert.ok(volcanoBox)
        assert.ok(Math.abs(volcanoBox.width - volcanoBox.height) < 2)
        assert.ok(volcanoBox.height <= page.viewportSize().height * 0.7)
        assert.equal(
          await page.locator('.de-statistics dt').first().evaluate((node) => getComputedStyle(node).overflowWrap),
          'normal',
        )
        await visible(page.getByText('123.40 CPM', { exact: false }))
        assert.equal(await page.getByText('Mean normalized count', { exact: false }).count(), 0)
        await visible(page.getByRole('link', { name: 'Download DEG results', exact: true }))
        assert.equal(await page.locator('.volcano-explorer-gene-set').evaluate((node) => node.open), true)
        await page.getByLabel('Highlight genes').fill('GENE3')
        assert.equal(otherPaths.some(wasRequested), false, 'Unselected model was prefetched')
        assert.ok(requests.some((path) => path.includes('/data/de-models/') && path.endsWith('.gz.bin')))
        await chooseMenu(page, 'Analysis', 'study_a · Age (per decade)')
        await page.locator('.volcano-explorer-about > summary').click()
        await visible(page.getByText('Effect: age · per year.', { exact: false }))
        assert.equal(await page.getByRole('link', { name: 'Download DEG results', exact: true }).count(), 0)
        assert.equal(wasRequested(other.recordsAsset), true, 'Selected model was not fetched')
        await visible(page.getByText('3 of 3 studies significant', { exact: true }))
        assert.ok((await page.locator('svg[aria-label^="Volcano plot"] title').allTextContents())
          .some((title) => title.includes('3 of 3 studies significant')))
      } finally {
        await context.close()
      }
    })

    await runJourney('release', 'shows static previews while DE and trajectory data load', async () => {
      const context = await browser.newContext()
      let resumeDe, resumeTrajectory
      const deGate = new Promise((resolve) => { resumeDe = resolve })
      const trajectoryGate = new Promise((resolve) => { resumeTrajectory = resolve })
      try {
        const page = await context.newPage()
        await page.route('**/data/de-models/**', async (route) => {
          await deGate
          await route.continue()
        })
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(page.locator('.de-preview img'))
        assert.match(await page.locator('.de-preview img').getAttribute('src'), /de-previews\/cd14-monocyte\.svg/)
        resumeDe()
        await visible(page.locator('.de-statistics'))

        await page.route('**/data/pipeline-trajectory-map.json*', async (route) => {
          await trajectoryGate
          await route.continue()
        })
        await choosePageFromGroup(page, 'More age analysis', 'Age-bin trajectories')
        await visible(page.locator('.trajectory-map-loading img'))
        resumeTrajectory()
        await visible(page.locator('.umap'))
        await visible(page.getByText('Each point is a cell type × gene pair. We ran UMAP on distances between their age-bin trajectory shapes.'))
        assert.equal(await page.locator('.trajectory-notice').count(), 0)
      } finally {
        resumeDe()
        resumeTrajectory()
        await context.close()
      }
    })

    await runJourney('release', 'uses legacy trajectory data when the compact map asset is absent', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        const legacy = {
          ...built,
          assets: built.assets.filter((asset) => asset.kind !== 'trajectoryMap' && !asset.path.startsWith('trajectory-details/')),
        }
        await page.route('**/data/release.json.gz.bin', (route) =>
          route.fulfill({ body: gzipSync(JSON.stringify(legacy)), contentType: 'application/octet-stream' }),
        )
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'More age analysis', 'Age-bin trajectories')
        await visible(page.locator('.umap'))
        assert.equal(await page.locator('.trajectory-notice').count(), 0)
        assert.equal(await page.getByText('Trajectory data is unavailable.').count(), 0)
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
          new URL(route.request().url()).pathname.endsWith('.gz.bin')
            ? route.fulfill({
                body: gzipSync(JSON.stringify(legacy)),
                contentType: 'application/octet-stream',
              })
            : route.fulfill({ json: legacy }),
        )
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        const plot = page.locator('svg[aria-label^="Volcano plot"]')
        await visible(plot)
        assert.equal(await plot.locator('circle').count(), 2)
        assert.equal(await plot.locator('title').filter({ hasText: 'GENE3' }).count(), 0)
        await page.getByLabel('Search genes').fill('GENE3')
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
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(page.locator('.de-statistics'))
        await page.getByLabel('Search genes').fill('GENE3')
        await chooseMenu(page, 'Matching genes (1)', 'GENE3')
        await visible(page.getByRole('heading', { name: 'GENE3', exact: true }))
        const statistics = page.locator('.de-statistics dd')
        await visible(statistics.nth(0).filter({ hasText: '2.30' }))
        await visible(statistics.nth(1).filter({ hasText: '0.50' }))
        await visible(statistics.nth(4).filter({ hasText: 'Unavailable' }))
        const plot = page.locator('svg[aria-label^="Volcano plot"]')
        await clickPoint(plot, 'GENE1: effect')
        await clickPoint(plot, 'GENE2: effect')
        await page.waitForFunction(() => {
          const selected = [...document.querySelectorAll('svg[aria-label^="Volcano plot"] circle')]
            .filter((circle) => circle.getAttribute('stroke') === '#b5453c')
          return selected.length === 1 && selected[0].querySelector('title')?.textContent?.startsWith('GENE2:')
        })
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'toggles volcano colors between expression and study support', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        const colorSwitch = page.getByLabel('Color points by reference expression')
        await visible(colorSwitch)
        assert.equal(await colorSwitch.isChecked(), false)
        await visible(page.getByText('Significant studies', { exact: true }))
        await visible(page.locator('.de-statistics').getByText('3 of 3 studies significant'))
        await visible(page.locator('.de-statistics').getByText('123.40 CPM'))
        assert.deepEqual(await page.locator('.de-support-endpoints span').allTextContents(), ['1', '–', '3'])
        await page.getByLabel('Highlight genes').fill('GENE1')
        await visible(page.locator('.de-volcano-plot circle[data-highlighted="true"]'))
        const pointOrder = await page.locator('.de-volcano-plot g').evaluateAll((groups) => groups
          .map((group) => {
            const point = group.querySelector('circle > title')?.textContent
            return point ? {
              highlighted: group.querySelector('circle[data-highlighted="true"]') !== null,
              selected: group.querySelector('circle[data-selected="true"]') !== null,
            } : null
          })
          .filter(Boolean))
        const highlightedIndices = pointOrder.flatMap((point, index) => point.highlighted ? [index] : [])
        const backgroundIndices = pointOrder.flatMap((point, index) => point.highlighted ? [] : [index])
        assert.ok(highlightedIndices.length > 0)
        assert.ok(Math.min(...highlightedIndices) > Math.max(...backgroundIndices), 'Gene-set points should be drawn above other points')
        assert.equal(pointOrder.at(-1).selected, true, 'The selected point should remain on top')
        assert.ok(await page.locator('.de-volcano-plot circle[aria-hidden="true"]').count() >= highlightedIndices.length)
        await colorSwitch.check()
        await visible(page.locator('.de-support-label').filter({ hasText: 'log₂(reference-bin CPM + 1)' }))
        await visible(page.locator('.de-statistics').getByText('3 of 3 studies significant'))
        await visible(page.locator('.de-statistics').getByText('123.40 CPM'))
        assert.deepEqual(await page.locator('.de-support-endpoints span').allTextContents(), ['5.85', '–', '6.96'])
        await assertLegendEndpointAlignment(page)
        await colorSwitch.uncheck()
        await visible(page.getByText('Significant studies', { exact: true }))
        await assertLegendEndpointAlignment(page)
        await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
        await chooseMenu(page, 'Analysis', 'Combined · Age (per decade)')
        await visible(page.locator('.de-volcano-plot circle').first())
        await visible(page.locator('.de-support-label').filter({ hasText: 'log₂(reference-bin CPM + 1)' }))
        assert.equal(await page.getByLabel('Color points by reference expression').isChecked(), true)
        assert.equal(await page.getByText('Significant studies', { exact: true }).count(), 0)
        assert.deepEqual(await page.locator('.de-support-endpoints span').allTextContents(), ['4.70', '–', '6.25'])
        assert.equal(await page.getByRole('heading', { name: 'Volcano plot', exact: true }).count(), 0)
        const pointFill = await page.locator('.de-volcano-plot circle').first().getAttribute('fill')
        assert.notEqual(pointFill, '#b9bdbc')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'retrieves non-significant retained trajectories for an explicit gene', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By gene')
        await visible(page.getByLabel('Gene symbol'))
        assert.equal(await page.locator('.trajectory-notice').count(), 0)
        await page.getByLabel('Gene symbol').fill('GENE2')
        await page
          .locator('.profile-scatter title')
          .filter({ hasText: 'Naive CD4 T cell' })
          .waitFor({ state: 'attached' })
        assert.equal(await page.locator('.profile-scatter circle').filter({ has: page.locator('title', { hasText: 'Naive CD4 T cell' }) }).getAttribute('fill'), '#597FC6')
        const profilePlot = page.locator('.profile-scatter')
        await clickPoint(profilePlot, 'CD14 monocyte')
        await clickPoint(profilePlot, 'Naive CD4 T cell')
        await page.waitForFunction(() => {
          const selected = [...document.querySelectorAll('.profile-scatter circle')]
            .filter((circle) => circle.getAttribute('stroke') === '#b5453c')
          return selected.length === 1 && selected[0].querySelector('title')?.textContent?.startsWith('Naive CD4 T cell')
        })
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'overview explains the workflow on desktop and mobile', async () => {
      const context = await browser.newContext({ viewport: { width: 375, height: 800 } })
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
        const workflow = page.getByRole('region', { name: 'How the data come together' })
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

    await runJourney('navigation', 'clock selector includes cross-cell-type mean predictions', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        const modelRequests = []
        page.on('request', (request) => {
          const url = new URL(request.url())
          if (url.pathname.includes('/data/clock-models/')) modelRequests.push(url.pathname)
        })
        await openApp(page, served.url)
        assert.equal(modelRequests.length, 0, 'Clock model records should not load before opening the clock page')
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
        await choosePageFromGroup(page, 'More age analysis', 'Transcriptomic clocks')
        await visible(page.getByRole('heading', { name: 'Predicting age from gene expression' }))
        await visible(page.getByLabel('Cell type: Naive CD4 T cell'))
        assert.equal(
          await page.getByLabel('Cell type: Naive CD4 T cell').locator('.choice-color').evaluate((dot) => getComputedStyle(dot).backgroundColor),
          'rgb(89, 127, 198)',
        )
        await visible(page.getByRole('heading', { name: 'Age-aware study-offset elastic net' }))
        await page.locator('.clock-method-card').first().waitFor({ state: 'visible' })
        assert.equal(await page.locator('.clock-method-card').count(), 1)
        assert.equal(await page.locator('.clock-plot svg').count(), 2)
        assert.equal(new Set(modelRequests).size, 1, 'Only the selected model should load its point data')
        assert.equal(await menuTrigger(page, 'Clock model').getAttribute('aria-label'), 'Clock model: Age-aware study-offset elastic net')
        await menuTrigger(page, 'Clock model').click()
        assert.deepEqual(
          await page.locator('.clock-model-select [role="menuitemradio"]').evaluateAll((buttons) =>
            buttons.map((button) => button.querySelector('span')?.textContent?.trim())),
          ['Elastic net', 'Equal-study elastic net', 'Rank elastic net', 'Study z-score elastic net', 'Age-aware study-offset elastic net'],
        )
        await page.keyboard.press('Escape')
        const viewportWidth = page.viewportSize().width
        const clockCard = await page.locator('.clock-method-card').first().boundingBox()
        const clockIntro = await page.locator('.clock-intro').boundingBox()
        assert.ok(clockCard && clockCard.width > viewportWidth * 0.75, 'Clock cards should use the wide centered content area')
        assert.ok(clockCard && Math.abs(clockCard.x + clockCard.width / 2 - viewportWidth / 2) < 2, 'Clock cards should be centered')
        assert.ok(clockIntro && clockIntro.width > 600, 'Clock introduction should not remain in a narrow text column')
        assert.equal(await page.locator('.clock-plot svg').first().evaluate((svg) => getComputedStyle(svg).backgroundColor), 'rgba(0, 0, 0, 0)')
        const axisTicks = await page.locator('.clock-plot svg').first().locator('.clock-tick').allTextContents()
        await chooseMenu(page, 'Clock model', 'Elastic net')
        await visible(page.getByRole('heading', { name: 'Elastic net', exact: true }))
        assert.deepEqual(await page.locator('.clock-plot svg').first().locator('.clock-tick').allTextContents(), axisTicks)
        assert.equal(new Set(modelRequests).size, 2, 'Changing models should fetch only the newly selected model')
        await chooseMenu(page, 'Cell type', 'Across cell types · mean clock')
        await visible(page.getByText('Cross-cell-type ensemble', { exact: true }))
        await visible(page.getByRole('heading', { name: 'Elastic net · mean across cell types' }))
        assert.equal(await page.locator('.clock-method-card').count(), 1)
        assert.equal(await page.locator('.clock-plot svg').count(), 2)
        await chooseMenu(page, 'Cell type', 'CD14 monocyte')
        await visible(page.getByRole('heading', { name: 'Age-aware study-offset elastic net', exact: true }))
        assert.equal(await page.locator('.clock-method-card').count(), 1)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        assert.equal(await menuTrigger(page, 'Cell type').getAttribute('aria-label'), 'Cell type: CD14 monocyte')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'single-cell landscape shows embedding and ranked markers', async () => {
      const context = await browser.newContext({ viewport: { width: 390, height: 844 } })
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
        await page.getByRole('button', { name: 'Cell states', exact: true }).click()
        await visible(page.getByRole('heading', { name: 'Cell states' }))
        await visible(page.getByRole('img', { name: 'UMAP of CD14 monocyte, colored by local cluster' }))
        assert.equal(await page.getByText('Explore the local embedding and clusters within each annotated cell type.').count(), 0)
        assert.equal(await page.getByText('Each point is a cell; colors and labels show local Leiden clusters.').count(), 0)
        assert.deepEqual(
          await page.locator('.single-cell-layout .single-cell-panel h2').evaluateAll((headings) =>
            headings.map((heading) => heading.textContent.trim()),
          ),
          [
            'Clusters in context',
            'Genes distinguishing clusters',
            'Age-DE genes on UMAP',
            'Cluster 0: fraction across age',
          ],
        )
        await visible(page.getByText('TESTMARKER', { exact: true }))
        await visible(page.getByText('OTHER', { exact: true }))
        await visible(page.getByRole('img', { name: /Cluster 0 sample-level fractions across age/ }))
        await visible(page.getByRole('cell', { name: '25.0%' }))
        await visible(page.getByRole('cell', { name: '+1.50' }))
        await visible(page.getByText('log₂ FC').first())
        await visible(page.getByText(/Dashed lines show study trends/))
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
        await openApp(page, served.url)
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
        await openApp(page, served.url)
        await page.evaluate(() => window.scrollTo(0, document.documentElement.scrollHeight))
        await page.locator('.entry-grid a[href="#composition"]').click()
        await visible(page.locator('.pipeline-composition'))
        await page.waitForFunction(() => window.scrollY === 0)
        assert.equal(new URL(page.url()).hash, '#composition')
        await page.getByRole('button', { name: 'Overview', exact: true }).click()
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
        await openApp(page, served.url)
        await visible(page.getByText('SYNTHETIC TEST DATA · EXPORT', { exact: false }))
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(page.getByText('Synthetic fixture: not biological evidence.'))
        await chooseMenu(page, 'Analysis', 'study_a · Age (per decade)')
        await page.locator('.volcano-explorer-about > summary').click()
        await visible(page.getByText('Effect: age · per year.', { exact: false }))
        await page.getByLabel('Search genes').fill('GENE2')
        await chooseMenu(page, 'Matching genes (1)', 'GENE2')
        await visible(page.getByRole('heading', { name: 'GENE2', exact: true }))
        await visible(page.locator('.de-statistics dd').first().filter({ hasText: '-2.10' }))
        await choosePageFromGroup(page, 'Differential expression', 'By gene')
        await visible(page.locator('.profiles-page'))
        assert.equal(await page.getByLabel('Gene symbol').inputValue(), 'GENE2')
        await visible(page.getByText('metadata unavailable', { exact: false }).first())
        await choosePageFromGroup(page, 'More age analysis', 'Age-bin trajectories')
        await visible(page.locator('.pair-card').getByRole('button', { name: 'GENE2', exact: true }))
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(page.getByText('Synthetic fixture: not biological evidence.'))
        await visible(page.locator('.fraction-chart-panel').first())
        const agePlot = page.locator('svg[aria-label^="Model-adjusted absolute fractions"]')
        let tooltip = await agePlot.locator('title').allTextContents()
        assert.ok(tooltip.some((value) => value.includes('25.00% estimated fraction')))
        await chooseMenu(page, 'Age plot denominator', 'Within retained PBMCs')
        await agePlot.locator('title').filter({ hasText: '40.00% estimated fraction' }).first().waitFor({ state: 'attached' })
        tooltip = await agePlot.locator('title').allTextContents()
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

    await runJourney('navigation', 'shows study publication details and analysis reports', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
        await page.getByRole('button', { name: 'Datasets', exact: true }).click()
        await visible(page.getByRole('heading', { name: 'study_a', exact: true }))
        assert.equal(
          await page.getByRole('link', { name: '10.1234/fixture' }).getAttribute('href'),
          'https://doi.org/10.1234/fixture',
        )
        await page.getByRole('button', { name: 'Modelling', exact: true }).click()
        await visible(page.getByRole('heading', { name: 'Models and interpretation' }))
        await visible(page.getByRole('heading', { name: 'Detailed model diagnostics' }))
        await visible(page.getByRole('link', { name: 'Open the report index' }))
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'supports keyboard selection and dismissal of menus', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
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
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(menuTrigger(page, 'Cell type'))
        const defaultCellType = (await menuTrigger(page, 'Cell type').getAttribute('aria-label')).replace(
          'Cell type: ',
          '',
        )
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(menuTrigger(page, 'Cell type'))
        await assertMenu(page, 'Cell type', defaultCellType)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(menuTrigger(page, 'Cell type'))
        await page.getByLabel('Search genes').fill('GENE2')
        await chooseMenu(page, 'Matching genes (1)', 'GENE2')
        await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
        assert.equal(await menuTrigger(page, 'Cell type').locator('.choice-color').evaluate((dot) => dot.style.backgroundColor), 'rgb(89, 127, 198)')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(page.locator('.pipeline-composition'))
        await assertMenu(page, 'Cell type', 'Naive CD4 T cell')
        await chooseMenu(page, 'Cell type', 'CD14 monocyte')
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(page.getByRole('heading', { name: 'GENE2', exact: true }))
        await assertMenu(page, 'Cell type', 'CD14 monocyte')
        await chooseMenu(page, 'Cell type', 'CD16 monocyte')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(
          page.getByText('Fraction summaries are unavailable for CD16 monocyte. Showing CD14 monocyte.'),
        )
        await assertMenu(page, 'Cell type', 'CD14 monocyte')
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(menuTrigger(page, 'Cell type'))
        await assertMenu(page, 'Cell type', 'CD16 monocyte')
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await page.getByRole('button', { name: 'Use CD14 monocyte' }).click()
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
        await visible(menuTrigger(page, 'Cell type'))
        await assertMenu(page, 'Cell type', 'CD14 monocyte')
      } finally {
        await context.close()
      }
    })

    await runJourney('navigation', 'shares the cell type selected from the trajectory map', async () => {
      const context = await browser.newContext()
      try {
        const page = await context.newPage()
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'More age analysis', 'Age-bin trajectories')
        const card = page.locator('.pair-card')
        await visible(card)
        const cellType = (await card.locator('.pair-card-cell-link').textContent()).trim()
        await card.locator('.pair-card-cell-link').click()
        await visible(page.locator('.de-page'))
        await page.getByRole('button', { name: 'Composition', exact: true }).click()
        await visible(page.locator('.pipeline-composition'))
        await assertMenu(page, 'Cell type', cellType)
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
          await openApp(page, served.url)
          await choosePageFromGroup(page, 'Differential expression', 'By cell type')
          await visible(menuTrigger(page, 'Analysis'))
          await chooseMenu(page, 'Analysis', 'study_a · Age (per decade)')
          await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
          await assertMenu(page, 'Analysis', 'study_a · Age (per decade)')
          await page.getByLabel('Search genes').fill('GENE2')
          await chooseMenu(page, 'Matching genes (1)', 'GENE2')
          await visible(page.locator('.de-statistics dd').first().filter({ hasText: '-0.67' }))
          await chooseMenu(page, 'Cell type', 'CD16 monocyte')
          await visible(
            page.getByText('The selected model is unavailable for CD16 monocyte.', { exact: false }),
          )
          await assertMenu(page, 'Analysis', 'study_a · Age (per decade) (unavailable)')
          await chooseMenu(page, 'Cell type', 'CD14 monocyte')
          await visible(page.locator('.de-statistics dd').first().filter({ hasText: '-2.10' }))
          await page.getByRole('button', { name: 'Composition', exact: true }).click()
          await visible(menuTrigger(page, 'Age plot denominator'))
          await chooseMenu(page, 'Age plot denominator', 'Within retained PBMCs')
          await chooseMenu(page, 'Cell type', 'Naive CD4 T cell')
          await assertMenu(page, 'Age plot denominator', 'Within retained PBMCs')
          await choosePageFromGroup(page, 'Differential expression', 'By cell type')
          await visible(menuTrigger(page, 'Analysis'))
          await assertMenu(page, 'Analysis', 'study_a · Age (per decade)')
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
        await openApp(page, served.url)
        await choosePageFromGroup(page, 'Differential expression', 'By cell type')
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
            new URL(route.request().url()).pathname.endsWith('.gz.bin')
              ? route.fulfill({
                  body: gzipSync(JSON.stringify(updated)),
                  contentType: 'application/octet-stream',
                })
              : route.fulfill({ json: updated }),
          )
          await page.route('**/data/gene-metadata.json*', (route) =>
            route.fulfill({ json: { releaseId: built.releaseId, schemaVersion: 1, genes: 'invalid' } }),
          )
          await openApp(page, served.url)
          await choosePageFromGroup(page, 'Differential expression', 'By gene')
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
            new URL(route.request().url()).pathname.endsWith('.gz.bin')
              ? route.fulfill({
                  body: gzipSync(JSON.stringify(payload)),
                  contentType: 'application/octet-stream',
                })
              : route.fulfill({ json: payload }),
          )
          await openApp(page, served.url)
          await choosePageFromGroup(page, 'Differential expression', 'By cell type')
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
