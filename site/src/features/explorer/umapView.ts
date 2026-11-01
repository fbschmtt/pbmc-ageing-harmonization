import type { OrthographicViewState } from '@deck.gl/core'
import type { TrajectoryRecord } from '../../data/atlas'

export type UmapViewState = Pick<OrthographicViewState, 'target' | 'zoom'>

export const entryZoomIncrease = 2

export function initialUmapView(
  fittedView: OrthographicViewState,
  selected: Pick<TrajectoryRecord, 'x' | 'y'> | null,
  savedView: UmapViewState | null,
): OrthographicViewState {
  if (savedView) return { ...fittedView, ...savedView }
  if (!selected) return fittedView
  return {
    ...fittedView,
    target: [selected.x, selected.y, 0],
    zoom: Number(fittedView.zoom ?? 0) + entryZoomIncrease,
  }
}
