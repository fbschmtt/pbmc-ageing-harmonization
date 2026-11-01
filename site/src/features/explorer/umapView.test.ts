import type { OrthographicViewState } from '@deck.gl/core'
import { describe, expect, it } from 'vitest'
import { entryZoomIncrease, initialUmapView } from './umapView'

const fittedView = {
  target: [0, 0, 0],
  zoom: 3,
} as OrthographicViewState

describe('initialUmapView', () => {
  it('centers a card-selected point and increases the fitted zoom', () => {
    const view = initialUmapView(fittedView, { x: 12, y: -4 }, null)
    expect(view.target).toEqual([12, -4, 0])
    expect(view.zoom).toBe(3 + entryZoomIncrease)
  })

  it('uses a saved view instead of recentering the selected point', () => {
    const view = initialUmapView(fittedView, { x: 12, y: -4 }, { target: [8, 9, 0], zoom: 6 })
    expect(view.target).toEqual([8, 9, 0])
    expect(view.zoom).toBe(6)
  })

  it('uses the global fitted view without a selection or saved view', () => {
    expect(initialUmapView(fittedView, null, null)).toEqual(fittedView)
  })
})
