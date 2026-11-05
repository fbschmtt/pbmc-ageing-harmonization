import { test } from 'node:test'
import assert from 'node:assert/strict'
import { mkdtemp, mkdir, readFile, rm, writeFile } from 'node:fs/promises'
import { tmpdir } from 'node:os'
import { join } from 'node:path'
import { createHash } from 'node:crypto'
import { createRelease } from './release.mjs'
import { verifyRelease } from './verify-release.mjs'
import { parseTrajectoryMetadata } from './trajectory-metadata.mjs'

const sha = (text) => createHash('sha256').update(text).digest('hex')
test('legacy trajectory metadata handles non-finite literals without changing strings', () => {
  assert.deepEqual(
    parseTrajectoryMetadata('{"text":"NaN, \\"Infinity\\"", "values":[NaN, Infinity, -Infinity]}'),
    { text: 'NaN, "Infinity"', values: [null, null, null] },
  )
  assert.throws(() => parseTrajectoryMetadata('{"value": NaNonsense}'), SyntaxError)
})
test('release validates and inventories a snapshot; failure preserves previous release', async () => {
  const root = await mkdtemp(join(tmpdir(), 'atlas-release-'))
  const previousExporterRevision = process.env.SITE_EXPORTER_REVISION
  process.env.SITE_EXPORTER_REVISION = 'contract-test-exporter-revision'
  try {
    const input = join(root, 'analysis/per_type'),
      type = join(input, 'monocyte/differential_expression')
    const target = join(root, 'data')
    await mkdir(type, { recursive: true })
    const options = { pathwayDirectory: join(root, 'pathways'), pantherPath: join(root, 'absent-panther') }
    await mkdir(options.pathwayDirectory)
    await writeFile(
      join(type, 'run_metadata.json'),
      JSON.stringify({
        cell_type: 'Monocyte',
        analysis_mode: 'test_only',
        interpretation_warning: 'Synthetic data',
        models: [
          {
            model: 'merged',
            status: 'complete',
            design: '~ age',
            n_samples: 10,
            study_sample_counts: { A: 10 },
            results_by_covariate: { age: 'merged.csv' },
          },
          { model: 'per_study', status: 'skipped', reason: 'insufficient samples' },
        ],
      }),
    )
    await writeFile(
      join(type, 'merged.csv'),
      'gene,log2FoldChange,padj,pvalue,lfcSE,baseMean,contrast\nGENE1,0.1,0.01,0.001,0.02,5,per year\nGENE2,0.2,0.5,0.4,0.03,7,per year\n',
    )
    await writeFile(join(type, 'report.html'), '<html><body>DE diagnostics</body></html>')
    const release = await createRelease(input, target, options)
    assert.equal(release.exporterRevision, 'contract-test-exporter-revision')
    assert.equal(release.inventory[0].models[1].status, 'skipped')
    const report = release.assets.find((item) => item.kind === 'report')
    assert.ok(report.path.startsWith(`reports/${release.releaseId}/`))
    assert.equal(
      await readFile(join(target, report.path), 'utf8'),
      '<html><body>DE diagnostics</body></html>',
    )
    assert.equal(report.warning, 'Synthetic data')
    assert.equal(release.inventory[0].models[1].reason, 'insufficient samples')
    const de = await readFile(join(target, 'de-index.json'), 'utf8')
    assert.equal(JSON.parse(de).releaseId, release.releaseId)
    assert.equal(release.assets.find((a) => a.path === 'de-index.json').sha256, sha(de))
    assert.equal(JSON.parse(de).cellTypes[0].warning, 'Synthetic data')
    const dePreview = JSON.parse(de).cellTypes[0].previewAsset
    assert.equal(dePreview, 'de-previews/monocyte.svg')
    assert.match(await readFile(join(target, dePreview), 'utf8'), /<svg xmlns=/)
    assert.ok(release.assets.some((asset) => asset.path === dePreview && asset.kind === 'de_preview'))
    const model = JSON.parse(de).cellTypes[0].models[0]
    assert.deepEqual(model.records, [])
    assert.equal(model.downloadAsset, 'de-downloads/monocyte/merged-age.csv')
    assert.equal(release.assets.find((asset) => asset.path === model.downloadAsset).kind, 'de_download')
    assert.match(await readFile(join(target, model.downloadAsset), 'utf8'), /GENE2,0\.2,0\.5/)
    const chunk = JSON.parse(await readFile(join(target, model.recordsAsset), 'utf8'))
    assert.deepEqual(
      chunk.cellTypes[0].models[0].records.map((record) => record.gene),
      ['GENE1'],
    )
    const extra = JSON.parse(await readFile(join(target, Object.values(model.recordBuckets)[0]), 'utf8'))
    assert.deepEqual(
      extra.cellTypes[0].models[0].records.map((record) => record.gene),
      ['GENE2'],
    )
    await verifyRelease(target)
    await writeFile(join(target, 'de-index.json'), de + ' ')
    await assert.rejects(verifyRelease(target), /checksum mismatch/)
    await writeFile(join(target, 'de-index.json'), de)
    const oldManifest = await readFile(join(target, 'release.json'), 'utf8')
    await writeFile(join(type, 'merged.csv'), 'gene,padj\nGENE1,0.01\n')
    await assert.rejects(createRelease(input, target, options), /missing log2FoldChange/)
    assert.equal(await readFile(join(target, 'release.json'), 'utf8'), oldManifest)
    assert.equal(await readFile(join(target, 'de-index.json'), 'utf8'), de)
    await writeFile(join(type, 'merged.csv'), 'gene,log2FoldChange,padj\nGENE1,0.1,0.01\n')
    const trajectory = join(root, 'analysis/combined/trajectory')
    await mkdir(trajectory, { recursive: true })
    await writeFile(join(trajectory, 'analysis_metadata.json'), '{}')
    await assert.rejects(createRelease(input, target, options), /Trajectory inputs are incomplete/)
    assert.equal(await readFile(join(target, 'release.json'), 'utf8'), oldManifest)
    await writeFile(join(trajectory, 'cross_cell_type_trajectory_clusters.csv'), 'gene\n')
    await writeFile(
      join(trajectory, 'analysis_metadata.json'),
      '{"common_bins":[],"cell_type_support":[{"cell_type":"Monocyte","age_bin_counts":{"20-30":NaN}}],"clustering":{"noise_fraction":NaN}}',
    )
    const legacy = await createRelease(input, target, options)
    assert.equal(legacy.context.ageSupport.cellTypes[0].ageBinCounts['20-30'], null)
    await verifyRelease(target)
  } finally {
    if (previousExporterRevision === undefined) delete process.env.SITE_EXPORTER_REVISION
    else process.env.SITE_EXPORTER_REVISION = previousExporterRevision
    await rm(root, { recursive: true, force: true })
  }
})
