import { describe, expect, it, vi } from 'vitest'
import { installOfflineRoutes } from './network.mjs'

describe('offline browser routes', () => {
  it('falls through for local requests so journey-specific routes can handle them', async () => {
    let handler
    const page = { route: vi.fn((_pattern, callback) => { handler = callback }) }
    await installOfflineRoutes(page, 'http://127.0.0.1:8000/')
    const fallback = vi.fn()
    await handler({ request: () => ({ url: () => 'http://127.0.0.1:8000/data/de-index.json' }), fallback })
    expect(fallback).toHaveBeenCalledOnce()
  })

  it('mocks Google Fonts and blocks other external requests', async () => {
    let handler
    const page = { route: vi.fn((_pattern, callback) => { handler = callback }) }
    await installOfflineRoutes(page, 'http://127.0.0.1:8000/')
    const font = { request: () => ({ url: () => 'https://fonts.googleapis.com/css2?family=Manrope' }), fulfill: vi.fn() }
    await handler(font)
    expect(font.fulfill).toHaveBeenCalledWith({ status: 200, contentType: 'text/css', body: '' })
    const external = { request: () => ({ url: () => 'https://example.com/data.json' }), abort: vi.fn() }
    await handler(external)
    expect(external.abort).toHaveBeenCalledOnce()
  })
})
