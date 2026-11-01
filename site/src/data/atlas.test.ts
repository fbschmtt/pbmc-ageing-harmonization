import { afterEach, beforeEach, describe, expect, it, vi } from 'vitest'

const payload = {
  schemaVersion: 1,
  releaseId: 'release-a',
  release: 'test',
  records: [],
  tissues: [],
  cellTypes: [],
  genes: [],
}
const manifest = {
  generatedAt: '2026-01-01',
  schemaVersion: 1,
  releaseId: 'release-a',
  assets: [{ path: 'pipeline-trajectory-index.json', kind: 'trajectory', bytes: 1, sha256: 'x' }],
}
beforeEach(() => vi.stubGlobal('DecompressionStream', undefined))
afterEach(() => {
  vi.unstubAllGlobals()
  vi.resetModules()
})
describe('trajectory data loading', () => {
  it('loads the pipeline release', async () => {
    const fetch = vi
      .fn()
      .mockResolvedValueOnce({ ok: true, json: async () => manifest })
      .mockResolvedValueOnce({ ok: true, json: async () => payload })
    vi.stubGlobal('fetch', fetch)
    const { loadAtlas } = await import('./atlas')
    expect(await loadAtlas()).toEqual(payload)
    expect(fetch.mock.calls[1][0]).toContain('pipeline-trajectory-index.json')
  })
  it('coalesces concurrent loads and retries failed assets', async () => {
    const fetch = vi
      .fn()
      .mockResolvedValueOnce({ ok: true, json: async () => manifest })
      .mockResolvedValueOnce({ ok: false })
      .mockResolvedValueOnce({ ok: true, json: async () => payload })
    vi.stubGlobal('fetch', fetch)
    const { loadAtlas } = await import('./atlas')
    const first = await Promise.allSettled([loadAtlas(), loadAtlas()])
    expect(first.map((result) => result.status)).toEqual(['rejected', 'rejected'])
    expect(fetch).toHaveBeenCalledTimes(2)
    expect(await loadAtlas()).toEqual(payload)
    expect(await loadAtlas()).toEqual(payload)
    expect(fetch).toHaveBeenCalledTimes(3)
  })
  it('rejects assets from a different release', async () => {
    vi.stubGlobal(
      'fetch',
      vi
        .fn()
        .mockResolvedValueOnce({ ok: true, json: async () => manifest })
        .mockResolvedValueOnce({ ok: true, json: async () => ({ ...payload, releaseId: 'old' }) }),
    )
    const { loadAtlas } = await import('./atlas')
    await expect(loadAtlas()).rejects.toThrow('does not match')
  })
  it('reports a missing manifest', async () => {
    vi.stubGlobal('fetch', vi.fn().mockResolvedValue({ ok: false }))
    const { loadAtlas } = await import('./atlas')
    await expect(loadAtlas()).rejects.toThrow('unavailable')
  })
})
