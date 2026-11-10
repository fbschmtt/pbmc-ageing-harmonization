import { gzipSync } from 'node:zlib'
import { afterEach, expect, it, vi } from 'vitest'
import { loadDifferentialDownload } from './differentialDownload'

const { manifest } = vi.hoisted(() => ({ manifest: vi.fn() }))
vi.mock('./release', async (original) => ({
  ...await original<typeof import('./release')>(),
  loadReleaseManifest: manifest,
}))
afterEach(() => { vi.unstubAllGlobals(); vi.resetAllMocks() })

function setup(path = 'de-downloads/cell/merged-age.csv.gz') {
  manifest.mockResolvedValue({ releaseId: 'release-1', assets: [{ path, kind: 'de_download' }] })
  const signal = new AbortController().signal
  return { path, signal }
}

it('extracts the original CSV bytes and metadata without changing precision or missing values', async () => {
  const { path, signal } = setup()
  const csv = '# Cell type: T célł\ngene,log2FoldChange,padj\nGENE1,0.12345678901234567,0.001\nGENE2,,\n'
  const fetcher = vi.fn().mockResolvedValue(new Response(new Uint8Array(gzipSync(csv))))
  vi.stubGlobal('fetch', fetcher)
  const blob = await loadDifferentialDownload(path, 'release-1', signal)
  expect(await blob.text()).toBe(csv)
  expect(blob.type).toBe('text/csv;charset=utf-8')
  expect(fetcher).toHaveBeenCalledWith(`/data/${path}?release=release-1`, { signal })
})

it('continues to support plain CSV downloads from older releases', async () => {
  const { path, signal } = setup('de-downloads/cell/merged-age.csv')
  vi.stubGlobal('DecompressionStream', undefined)
  vi.stubGlobal('fetch', vi.fn().mockResolvedValue(new Response('gene,padj\nGENE1,0.01\n')))
  expect(await (await loadDifferentialDownload(path, 'release-1', signal)).text()).toBe('gene,padj\nGENE1,0.01\n')
})

it('rejects missing downloads, HTTP failures and corrupt or truncated gzip data', async () => {
  const { path, signal } = setup()
  const fetcher = vi.fn()
  vi.stubGlobal('fetch', fetcher)
  await expect(loadDifferentialDownload('missing.csv.gz', 'release-1', signal)).rejects.toThrow('not included')
  expect(fetcher).not.toHaveBeenCalled()
  fetcher.mockResolvedValueOnce(new Response(null, { status: 404 }))
  await expect(loadDifferentialDownload(path, 'release-1', signal)).rejects.toThrow('unavailable')
  for (const bytes of [new TextEncoder().encode('invalid gzip'), new Uint8Array(gzipSync('gene,padj\n').subarray(0, -5))]) {
    fetcher.mockResolvedValueOnce(new Response(bytes))
    await expect(loadDifferentialDownload(path, 'release-1', signal)).rejects.toThrow()
  }
})

it('rejects downloads belonging to a different release before requesting bytes', async () => {
  const { path, signal } = setup()
  const fetcher = vi.fn()
  vi.stubGlobal('fetch', fetcher)
  await expect(loadDifferentialDownload(path, 'old-release', signal)).rejects.toThrow('different release')
  expect(fetcher).not.toHaveBeenCalled()
})
