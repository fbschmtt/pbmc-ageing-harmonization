import { useCallback, useEffect, useMemo, useRef, useState, type Dispatch } from 'react'
import { loadAtlas, loadTrajectoryChunk, loadTrajectoryMap, type Atlas, type TrajectoryMap, type TrajectoryPoint, type TrajectoryRecord } from '../data/atlas'
import { geneBucket } from '../data/geneBucket.mjs'
import { parseGeneSet } from '../data/geneSet'
import { loadGeneMetadata, type GeneMetadataIndex } from '../data/geneMetadata'
import { loadPathwaySets, type PathwaySet } from '../data/pathways'
import {
  assetFailureKind,
  loadReleaseManifest,
  type ReleaseManifest,
  type AssetFailureKind,
} from '../data/release'
import type { Page } from './navigation'
import type { AtlasSessionAction, AtlasSessionState } from './session'

export function useSiteData(page: Page, dispatch: Dispatch<AtlasSessionAction>, session: AtlasSessionState) {
  const [releaseLoading, setReleaseLoading] = useState(true)
  const [release, setRelease] = useState<ReleaseManifest | null>(null)
  const [atlas, setAtlas] = useState<Atlas | null>(null)
  const [trajectoryMap, setTrajectoryMap] = useState<TrajectoryMap | null>(null)
  const [loadError, setLoadError] = useState<AssetFailureKind | null>(null)
  const [annotationError, setAnnotationError] = useState<AssetFailureKind | null>(null)
  const [pathwayError, setPathwayError] = useState<AssetFailureKind | null>(null)
  const [pathwaySets, setPathwaySets] = useState<PathwaySet[]>([])
  const [geneMetadata, setGeneMetadata] = useState<GeneMetadataIndex | null>(null)
  const [extraRecords, setExtraRecords] = useState<Record<string, TrajectoryRecord[]>>({})
  const pendingChunks = useRef(new Set<string>())
  const [extraError, setExtraError] = useState(false)
  const requestChunk = useCallback((path: string) => {
    if (extraRecords[path] || pendingChunks.current.has(path)) return
    pendingChunks.current.add(path)
    setExtraError(false)
    loadTrajectoryChunk(path)
      .then((chunk) => setExtraRecords((previous) => ({ ...previous, [path]: chunk.records })))
      .catch(() => setExtraError(true))
      .finally(() => pendingChunks.current.delete(path))
  }, [extraRecords])
  const ensureTrajectoryDetails = useCallback((point: TrajectoryPoint) => {
    if (atlas?.records.some((record) => record.id === point.id)) return
    const path = trajectoryMap?.geneBuckets?.[geneBucket(point.gene)]
    if (path) requestChunk(path)
  }, [atlas, trajectoryMap, requestChunk])
  const requestedText =
    page === 'geneTypeMatrix'
      ? session.highlightGenes
      : page === 'genes'
        ? session.geneProfile.gene || session.selectedGene || ''
        : page === 'cellTypes'
          ? session.selectedGene || ''
          : page === 'trajectories'
            ? session.trajectory.query || (session.trajectory.selectedId ? session.selectedGene || '' : '')
            : ''
  const requestedGenes = useMemo(() => parseGeneSet(requestedText), [requestedText])
  useEffect(() => {
    const source = page === 'trajectories' ? trajectoryMap : atlas
    if (!source?.geneBuckets) return
    const paths = [
      ...new Set(
        [...requestedGenes]
          .map((gene) => source.geneBuckets?.[geneBucket(gene)])
          .filter((path): path is string => !!path && !extraRecords[path]),
      ),
    ]
    if (!paths.length) return
    const timer = window.setTimeout(() => {
      for (const path of paths) requestChunk(path)
    }, 250)
    return () => {
      window.clearTimeout(timer)
    }
  }, [atlas, trajectoryMap, page, requestedGenes, extraRecords, requestChunk])
  const visibleAtlas = useMemo(() => {
    if (!atlas?.geneBuckets || !requestedGenes.size) return atlas
    const additional = Object.values(extraRecords)
      .flat()
      .filter((record) => requestedGenes.has(record.gene.toUpperCase()))
    return { ...atlas, records: [...atlas.records, ...additional] }
  }, [atlas, extraRecords, requestedGenes])
  const visibleTrajectoryMap = useMemo(() => {
    if (!trajectoryMap) return null
    const baseIds = new Set(trajectoryMap.records.map((record) => record.id))
    const additional = Object.values(extraRecords).flat()
      .filter((record) => requestedGenes.has(record.gene.toUpperCase()) && !baseIds.has(record.id))
      .map(({ id, x, y, cellType, tissue, gene, adjustedPValue, cluster }): TrajectoryPoint =>
        ({ id, x, y, cellType, tissue, gene, adjustedPValue, cluster }))
    return { ...trajectoryMap, records: [...trajectoryMap.records, ...additional] }
  }, [trajectoryMap, extraRecords, requestedGenes])
  useEffect(() => {
    const source = page === 'trajectories' ? trajectoryMap : atlas
    if (!source?.geneBuckets || !['genes', 'cellTypes', 'trajectories', 'geneTypeMatrix'].includes(page))
      return
    const genes = parseGeneSet(session.highlightGenes)
    if (!genes.size) return
    const paths = [
      ...new Set(
        [...genes]
          .map((gene) => source.geneBuckets?.[geneBucket(gene)])
          .filter((path): path is string => !!path),
      ),
    ]
    const timer = window.setTimeout(() => {
      for (const path of paths) requestChunk(path)
    }, 1000)
    return () => {
      window.clearTimeout(timer)
    }
  }, [atlas, trajectoryMap, session.highlightGenes, page, extraRecords, requestChunk])
  useEffect(() => {
    if (!['genes', 'cellTypes', 'trajectories'].includes(page) || geneMetadata) return
    let active = true
    loadGeneMetadata()
      .then((value) => {
        if (active) {
          setGeneMetadata(value)
          setAnnotationError(null)
        }
      })
      .catch((error) => {
        if (active) {
          setGeneMetadata(null)
          setAnnotationError(assetFailureKind(error))
        }
      })
    return () => {
      active = false
    }
  }, [geneMetadata, page])
  useEffect(() => {
    let active = true
    loadReleaseManifest()
      .then((value) => {
        if (active) {
          setRelease(value)
          setReleaseLoading(false)
        }
      })
      .catch(() => {
        if (active) {
          setRelease(null)
          setReleaseLoading(false)
        }
      })
    return () => {
      active = false
    }
  }, [])
  useEffect(() => {
    if (page !== 'trajectories' || trajectoryMap) return
    let active = true
    loadTrajectoryMap()
      .then((value) => {
        if (active) setTrajectoryMap(value)
      })
      .catch((error) => {
        if (!active) return
        if (assetFailureKind(error) !== 'unavailable') {
          setLoadError(assetFailureKind(error))
          return
        }
        loadAtlas()
          .then((legacy) => {
            if (!active) return
            setAtlas(legacy)
            setTrajectoryMap({
              schemaVersion: 1,
              release: legacy.release,
              warning: legacy.warning,
              description: legacy.description,
              records: legacy.records.map(({ id, x, y, cellType, tissue, gene, adjustedPValue, cluster }) =>
                ({ id, x, y, cellType, tissue, gene, adjustedPValue, cluster })),
              tissues: legacy.tissues,
              cellTypes: legacy.cellTypes,
              genes: legacy.genes,
              geneBuckets: legacy.geneBuckets,
            })
          })
          .catch((fallbackError) => {
            if (active) setLoadError(assetFailureKind(fallbackError))
          })
      })
    return () => {
      active = false
    }
  }, [page, trajectoryMap])
  useEffect(() => {
    if (!['genes', 'cellTypes', 'geneTypeMatrix'].includes(page) || atlas) return
    let active = true
    loadAtlas()
      .then((value) => {
        if (active) setAtlas(value)
      })
      .catch((error) => {
        if (active) setLoadError(assetFailureKind(error))
      })
    return () => {
      active = false
    }
  }, [page, atlas])
  useEffect(() => {
    if (!['pathways', 'differential', 'geneTypeMatrix', 'trajectories', 'cellTypes'].includes(page)) return
    let active = true
    loadPathwaySets()
      .then((sets) => {
        if (active) {
          setPathwaySets(sets)
          dispatch({ type: 'initializePathway', pathways: sets })
        }
      })
      .catch((error) => {
        if (active) {
          setPathwaySets([])
          setPathwayError(assetFailureKind(error))
        }
      })
    return () => {
      active = false
    }
  }, [dispatch, page])
  return {
    release,
    releaseLoading,
    atlas: visibleAtlas,
    trajectoryMap: visibleTrajectoryMap,
    trajectoryDetails: [...(atlas?.records ?? []), ...Object.values(extraRecords).flat()],
    ensureTrajectoryDetails,
    loadError,
    pathwaySets,
    geneMetadata,
    annotationError,
    pathwayError,
    extraError,
  }
}
