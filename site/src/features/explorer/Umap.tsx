import { OrthographicView, type OrthographicViewState, type PickingInfo } from '@deck.gl/core'
import { ScatterplotLayer } from '@deck.gl/layers'
import DeckGL from '@deck.gl/react'
import { useCallback, useEffect, useLayoutEffect, useMemo, useRef, useState } from 'react'
import type { TrajectoryPoint } from '../../data/atlas'
import { deckColors } from '../../theme/visualEncoding'
import { initialUmapView, type UmapViewState } from './umapView'

type UmapProps = {
  records: TrajectoryPoint[]
  matchingRecords: TrajectoryPoint[]
  selected: TrajectoryPoint | null
  focusSelection: boolean
  highlightedGenes: ReadonlySet<string>
  onSelect: (record: TrajectoryPoint) => void
  initialViewState: UmapViewState | null
  resetVersion: number
  onViewStateChange: (viewState: UmapViewState) => void
}

type Size = { width: number; height: number }

const view = new OrthographicView({ id: 'umap', flipY: false })
const getPosition = (record: TrajectoryPoint): [number, number] => [record.x, record.y]
const getFillColor = (record: TrajectoryPoint) =>
  record.tissue === 'brain' ? deckColors.brain : deckColors.blood
const getSubduedFillColor = (record: TrajectoryPoint) =>
  record.tissue === 'brain' ? deckColors.subduedBrain : deckColors.subduedBlood
const getHighlightedFillColor = (record: TrajectoryPoint) =>
  record.tissue === 'brain' ? deckColors.highlightedBrain : deckColors.highlightedBlood

function fitView(records: TrajectoryPoint[], { width, height }: Size): OrthographicViewState {
  if (records.length === 0) return { target: [0, 0, 0], zoom: 0 }

  let minX = Infinity
  let maxX = -Infinity
  let minY = Infinity
  let maxY = -Infinity
  for (const record of records) {
    minX = Math.min(minX, record.x)
    maxX = Math.max(maxX, record.x)
    minY = Math.min(minY, record.y)
    maxY = Math.max(maxY, record.y)
  }
  const xSpan = Math.max(maxX - minX, 1)
  const ySpan = Math.max(maxY - minY, 1)
  const padding = 48
  const scale = Math.min(Math.max(width - padding * 2, 1) / xSpan, Math.max(height - padding * 2, 1) / ySpan)

  return {
    target: [(minX + maxX) / 2, (minY + maxY) / 2, 0],
    zoom: Math.log2(Math.max(scale, Number.EPSILON)),
  }
}

export function Umap({
  records,
  matchingRecords,
  selected,
  focusSelection,
  highlightedGenes,
  onSelect,
  initialViewState,
  resetVersion,
  onViewStateChange,
}: UmapProps) {
  const containerRef = useRef<HTMLDivElement>(null)
  const onSelectRef = useRef(onSelect)
  const onViewStateChangeRef = useRef(onViewStateChange)
  const savedViewStateRef = useRef<UmapViewState | null>(initialViewState)
  const saveTimerRef = useRef<number | null>(null)
  const hasAppliedInitialFitRef = useRef(Boolean(initialViewState))
  const previousResetVersionRef = useRef(resetVersion)
  onSelectRef.current = onSelect
  onViewStateChangeRef.current = onViewStateChange
  const [size, setSize] = useState<Size>({ width: 640, height: 500 })
  const [renderError, setRenderError] = useState<string | null>(null)
  const fittedViewState = useMemo(() => fitView(records, size), [records, size])
  const entryViewState = useMemo<OrthographicViewState>(
    () => initialUmapView(fittedViewState, focusSelection ? selected : null, initialViewState),
    [fittedViewState, initialViewState, selected, focusSelection],
  )
  const [viewState, setViewState] = useState<OrthographicViewState>(() => entryViewState)

  const saveViewState = useCallback((nextViewState: OrthographicViewState, delay = 150) => {
    const savedViewState: UmapViewState = { target: nextViewState.target, zoom: nextViewState.zoom }
    savedViewStateRef.current = savedViewState
    if (saveTimerRef.current !== null) window.clearTimeout(saveTimerRef.current)
    saveTimerRef.current = window.setTimeout(() => {
      onViewStateChangeRef.current(savedViewState)
      saveTimerRef.current = null
    }, delay)
  }, [])

  useLayoutEffect(() => {
    const element = containerRef.current
    if (!element) return

    const updateSize = () => {
      const { width, height } = element.getBoundingClientRect()
      setSize({ width, height })
    }
    updateSize()
    const observer = new ResizeObserver(updateSize)
    observer.observe(element)
    return () => observer.disconnect()
  }, [])

  useEffect(() => {
    if (!hasAppliedInitialFitRef.current) {
      hasAppliedInitialFitRef.current = true
      setViewState(entryViewState)
      saveViewState(entryViewState, 0)
      return
    }
  }, [entryViewState, saveViewState])

  useEffect(() => {
    if (previousResetVersionRef.current === resetVersion) return
    previousResetVersionRef.current = resetVersion
    setViewState(fittedViewState)
    saveViewState(fittedViewState, 0)
  }, [fittedViewState, resetVersion, saveViewState])

  useEffect(
    () => () => {
      if (saveTimerRef.current !== null) window.clearTimeout(saveTimerRef.current)
      if (savedViewStateRef.current) onViewStateChangeRef.current(savedViewStateRef.current)
    },
    [],
  )

  // Keep points about 2.75 px at the fitted view, then let them grow with zoom.
  // A pixel maximum prevents extreme zoom levels from producing giant circles.
  const pointRadiusScale = 2.75 / 2 ** Number(fittedViewState.zoom ?? 0)
  const pointsLayer = useMemo(
    () =>
      new ScatterplotLayer<TrajectoryPoint>({
        id: 'umap-points',
        data: records,
        getPosition,
        getFillColor: getSubduedFillColor,
        getRadius: 1,
        radiusScale: pointRadiusScale,
        radiusUnits: 'common',
        radiusMinPixels: 2.75,
        radiusMaxPixels: 10,
        pickable: false,
      }),
    [pointRadiusScale, records],
  )

  const matchingPointsLayer = useMemo(
    () =>
      new ScatterplotLayer<TrajectoryPoint>({
        id: 'umap-matching-points',
        data: matchingRecords,
        getPosition,
        getFillColor,
        getRadius: 1,
        radiusScale: pointRadiusScale,
        radiusUnits: 'common',
        radiusMinPixels: 2.75,
        radiusMaxPixels: 10,
        pickable: true,
        autoHighlight: true,
        highlightColor: [29, 52, 48, 70],
        onClick: ({ object }: PickingInfo<TrajectoryPoint>) => {
          if (object) onSelectRef.current(object)
        },
      }),
    [matchingRecords, pointRadiusScale],
  )

  const geneHighlightLayer = useMemo(() => {
    const highlightedRecords = matchingRecords.filter((record) =>
      highlightedGenes.has(record.gene.toUpperCase()),
    )
    return highlightedRecords.length
      ? new ScatterplotLayer<TrajectoryPoint>({
          id: 'umap-gene-highlights',
          data: highlightedRecords,
          getPosition,
          getFillColor: getHighlightedFillColor,
          getLineColor: [0, 0, 0, 255],
          getRadius: 1,
          radiusScale: pointRadiusScale * 1.45,
          getLineWidth: 1.25,
          radiusUnits: 'common',
          radiusMinPixels: 4,
          radiusMaxPixels: 13,
          lineWidthUnits: 'pixels',
          stroked: true,
          pickable: true,
          onClick: ({ object }: PickingInfo<TrajectoryPoint>) => {
            if (object) onSelectRef.current(object)
          },
        })
      : null
  }, [highlightedGenes, matchingRecords, pointRadiusScale])

  const selectionLayer = useMemo(
    () =>
      selected
        ? new ScatterplotLayer<TrajectoryPoint>({
            id: 'umap-selection',
            data: [selected],
            getPosition,
            getFillColor,
            getLineColor: deckColors.selected,
            getRadius: 1,
            radiusScale: pointRadiusScale * 2,
            getLineWidth: 2,
            radiusUnits: 'common',
            radiusMinPixels: 6,
            radiusMaxPixels: 15,
            lineWidthUnits: 'pixels',
            stroked: true,
            pickable: false,
          })
        : null,
    [pointRadiusScale, selected],
  )

  const layers = useMemo(
    () => [pointsLayer, matchingPointsLayer, geneHighlightLayer, selectionLayer].filter(Boolean),
    [geneHighlightLayer, matchingPointsLayer, pointsLayer, selectionLayer],
  )

  if (records.length === 0) return <div className="empty-map">No trajectory records are available.</div>
  if (renderError)
    return (
      <div className="empty-map" role="alert">
        {renderError}
      </div>
    )

  return (
    <div
      ref={containerRef}
      className="umap umap-deck"
      role="img"
      aria-label={`Age-bin trajectories map with ${records.length.toLocaleString()} points, ${matchingRecords.length.toLocaleString()} matching the focus. Drag to pan and scroll to zoom.`}
    >
      <DeckGL
        views={view}
        viewState={viewState}
        controller={{ dragRotate: false }}
        layers={layers}
        onViewStateChange={({ viewState: nextViewState }) => {
          const next = nextViewState as OrthographicViewState
          setViewState(next)
          saveViewState(next)
        }}
        getCursor={({ isDragging, isHovering }) =>
          isDragging ? 'grabbing' : isHovering ? 'pointer' : 'grab'
        }
        getTooltip={({ object }: PickingInfo<TrajectoryPoint>) =>
          object
            ? {
                text: `${object.gene} · ${object.cellType}`,
                className: 'umap-tooltip',
              }
            : null
        }
        onError={() =>
          setRenderError(
            'The similarity map could not start accelerated rendering in this browser. Focus and gene-set controls remain available.',
          )
        }
      />
    </div>
  )
}
