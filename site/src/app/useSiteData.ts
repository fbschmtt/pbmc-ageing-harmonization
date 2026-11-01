import { useEffect, useMemo, useState, type Dispatch } from 'react'
import { loadAtlas, loadTrajectoryChunk, type Atlas, type TrajectoryRecord } from '../data/atlas'
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
  const [loadError, setLoadError] = useState<AssetFailureKind | null>(null)
  const [annotationError, setAnnotationError] = useState<AssetFailureKind | null>(null)
  const [pathwayError, setPathwayError] = useState<AssetFailureKind | null>(null)
  const [pathwaySets, setPathwaySets] = useState<PathwaySet[]>([])
  const [geneMetadata, setGeneMetadata] = useState<GeneMetadataIndex | null>(null)
  const [extraRecords, setExtraRecords] = useState<Record<string, TrajectoryRecord[]>>({})
  const [extraError, setExtraError] = useState(false)
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
    if (!atlas?.geneBuckets) return
    const paths = [
      ...new Set(
        [...requestedGenes]
          .map((gene) => atlas.geneBuckets?.[geneBucket(gene)])
          .filter((path): path is string => !!path && !extraRecords[path]),
      ),
    ]
    if (!paths.length) return
    let active = true
    const timer = window.setTimeout(() => {
      setExtraError(false)
      Promise.all(paths.map(async (path) => [path, (await loadTrajectoryChunk(path)).records] as const))
        .then((entries) => {
          if (active) setExtraRecords((previous) => ({ ...previous, ...Object.fromEntries(entries) }))
        })
        .catch(() => {
          if (active) setExtraError(true)
        })
    }, 250)
    return () => {
      active = false
      window.clearTimeout(timer)
    }
  }, [atlas, requestedGenes, extraRecords])
  const visibleAtlas = useMemo(() => {
    if (!atlas?.geneBuckets || !requestedGenes.size) return atlas
    const additional = Object.values(extraRecords)
      .flat()
      .filter((record) => requestedGenes.has(record.gene.toUpperCase()))
    return { ...atlas, records: [...atlas.records, ...additional] }
  }, [atlas, extraRecords, requestedGenes])
  useEffect(() => {
    if (!atlas?.geneBuckets || !['genes', 'cellTypes', 'trajectories', 'geneTypeMatrix'].includes(page))
      return
    const genes = parseGeneSet(session.highlightGenes)
    if (!genes.size) return
    const paths = [
      ...new Set(
        [...genes]
          .map((gene) => atlas.geneBuckets?.[geneBucket(gene)])
          .filter((path): path is string => !!path),
      ),
    ]
    let active = true
    let cursor = 0
    const worker = async () => {
      while (active && cursor < paths.length) {
        const path = paths[cursor++]
        try {
          await loadTrajectoryChunk(path)
        } catch {
          /* An explicit lookup reports errors. */
        }
      }
    }
    const timer = window.setTimeout(() => {
      void worker()
      void worker()
    }, 1000)
    return () => {
      active = false
      window.clearTimeout(timer)
    }
  }, [atlas, session.highlightGenes, page])
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
    if (!['genes', 'cellTypes', 'trajectories', 'geneTypeMatrix'].includes(page) || atlas) return
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
    loadError,
    pathwaySets,
    geneMetadata,
    annotationError,
    pathwayError,
    extraError,
  }
}
