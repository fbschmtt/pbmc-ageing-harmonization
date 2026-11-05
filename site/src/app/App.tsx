import { useEffect, useReducer } from 'react'
import type { PathwaySet } from '../data/pathways'
import { DifferentialPage } from '../features/differential/DifferentialPage'
import { Explorer } from '../features/explorer/Explorer'
import { PipelineComposition } from '../features/composition/PipelineComposition'
import { InformationalPage } from '../features/content/InformationalPage'
import { OverviewPage } from '../features/content/OverviewPage'
import { GeneTypeMatrix } from '../features/pathways/PathwaysPage'
import { PathwaySelector } from '../features/pathways/PathwaySelector'
import { ProfilesPage } from '../features/profiles/ProfilesPage'
import type { TrajectoryPoint, TrajectoryRecord } from '../data/atlas'
import { navigation } from './navigation'
import { usePageNavigation } from './usePageNavigation'
import { useSiteData } from './useSiteData'
import { releaseLabel } from './releaseLabel'
import { TrajectoryMenu } from './TrajectoryMenu'
import { atlasSessionReducer, initialAtlasSession } from './session'
import { SingleCellPage } from '../features/singleCell/SingleCellPage'
import { preloadCombinedAgeModel } from '../data/differential'

export function App() {
  const [session, dispatch] = useReducer(atlasSessionReducer, initialAtlasSession)
  useEffect(() => {
    const timer = window.setTimeout(() => {
      void preloadCombinedAgeModel(session.selectedCellType).catch(() => {})
    }, 1200)
    return () => window.clearTimeout(timer)
  }, [session.selectedCellType])
  const { page, navigate } = usePageNavigation()
  const {
    release,
    releaseLoading,
    atlas,
    trajectoryMap,
    trajectoryDetails,
    ensureTrajectoryDetails,
    loadError,
    pathwaySets,
    geneMetadata,
    annotationError,
    pathwayError,
    extraError,
  } = useSiteData(page, dispatch, session)
  const updateHighlightGenes = (value: string) => dispatch({ type: 'editGeneSet', value })
  const selectPathway = (pathway: PathwaySet) => dispatch({ type: 'selectPathway', pathway })
  const activePathway = pathwaySets.find((pathway) => pathway.id === session.activePathwayId) ?? null
  const dePreviewAssets = release?.assets.filter((asset) => asset.kind === 'de_preview') ?? []
  const selectedPreviewSlug = session.selectedCellType
    ?.toLowerCase()
    .replace(/[^a-z0-9]+/g, '-')
    .replace(/^-|-$/g, '')
  const loadingPreviewPath =
    dePreviewAssets.find((asset) => asset.path === `de-previews/${selectedPreviewSlug}.svg`)?.path ??
    dePreviewAssets.find((asset) => asset.path === 'de-previews/naive-cd8-t-cell.svg')?.path ??
    dePreviewAssets[0]?.path
  const trajectoryNotice = page === 'trajectories' ? trajectoryMap : atlas
  const openProfile = (record: TrajectoryRecord, tab: 'gene' | 'cellType', selectPair = true) => {
    dispatch({ type: 'selectGene', gene: record.gene, cellType: record.cellType })
    dispatch({ type: 'openProfile', tab, record, selectPair })
    navigate(tab === 'gene' ? 'genes' : 'cellTypes')
  }
  const openGene = (gene: string) => {
    const record = atlas?.records.find((item) => item.gene === gene)
    if (record) openProfile(record, 'gene', false)
  }
  const openCellType = (cellType: string) => {
    dispatch({ type: 'selectCellType', cellType })
    dispatch({
      type: 'setProfile',
      tab: 'cellType',
      state: { ...session.cellTypeProfile, cellType, selectedId: null },
    })
    navigate('cellTypes')
  }
  const openTrajectory = (record: TrajectoryRecord) => {
    dispatch({ type: 'selectGene', gene: record.gene, cellType: record.cellType })
    dispatch({ type: 'openTrajectory', record })
    navigate('trajectories')
  }
  const pairActions = {
    onOpenGene: (record: TrajectoryRecord) => openProfile(record, 'gene'),
    onOpenCellType: (record: TrajectoryRecord) => openProfile(record, 'cellType'),
    onOpenTrajectory: openTrajectory,
  }

  return (
    <main className={`atlas-main ${page === 'trajectories' ? 'is-trajectory-page' : ''}`}>
      <header className="site-header">
        <button className="brand" onClick={() => navigate('overview')}>
          <span className="brand-mark">A</span>
          <span>
            Ageing
            <br />
            Online Atlas
          </span>
        </button>
        <nav aria-label="Primary navigation">
          {navigation.map((item) => {
            if (item.id === 'genes') {
              return <TrajectoryMenu key="age-bin-trajectories" page={page} navigate={navigate} />
            }
            if (item.id === 'cellTypes' || item.id === 'trajectories' || item.id === 'geneTypeMatrix')
              return null
            return (
              <button
                key={item.id}
                onClick={() => navigate(item.id)}
                className={page === item.id ? 'active' : ''}
              >
                {item.label}
              </button>
            )
          })}
        </nav>
        <span className="release" title={release ? `Export ID: ${release.releaseId}` : undefined}>
          {releaseLoading ? 'LOADING RELEASE…' : releaseLabel(release)}
        </span>
      </header>
      {annotationError &&
        annotationError !== 'unavailable' &&
        ['genes', 'cellTypes', 'trajectories'].includes(page) && (
          <p role="status">
            Gene annotations could not be loaded because the release data is invalid or mismatched.
          </p>
        )}
      {extraError && <p role="status">Requested trajectory results could not be loaded.</p>}
      {pathwayError &&
        pathwayError !== 'unavailable' &&
        ['pathways', 'differential', 'geneTypeMatrix'].includes(page) && (
          <p role="status">
            Pathway data could not be loaded because the release data is invalid or mismatched.
          </p>
        )}
      {page === 'overview' && <OverviewPage onNavigate={navigate} />}
      {page === 'differential' && (
        <DifferentialPage
          loadingPreviewPath={loadingPreviewPath}
          releaseId={release?.releaseId}
          highlightGenes={session.highlightGenes}
          activePathway={activePathway}
          preferredGene={session.selectedGene}
          preferredCellType={session.selectedCellType}
          onSelectCellType={(cellType) => dispatch({ type: 'selectCellType', cellType })}
          preferredModelId={session.selectedDEModelId}
          onSelectModel={(modelId) => dispatch({ type: 'selectDEModel', modelId })}
          pathwaySets={pathwaySets}
          onSelectPathway={selectPathway}
          onHighlightGenesChange={updateHighlightGenes}
          onSelectGene={(gene, cellType) => dispatch({ type: 'selectGene', gene, cellType })}
        />
      )}
      {['genes', 'cellTypes', 'trajectories', 'geneTypeMatrix'].includes(page) && trajectoryNotice && (
        <div className="trajectory-notice">
          {trajectoryNotice.warning && <strong>{trajectoryNotice.warning} </strong>}
          {trajectoryNotice.description}
        </div>
      )}
      {!atlas && !trajectoryMap && ['genes', 'cellTypes', 'trajectories', 'geneTypeMatrix'].includes(page) && (
        <section className="content-page">
          <h1>
            {loadError === 'invalid'
              ? 'Trajectory data is invalid.'
              : loadError === 'mismatch'
                ? 'Trajectory data belongs to a different release.'
                : loadError
                  ? 'Trajectory data is unavailable.'
                  : 'Loading trajectory data…'}
          </h1>
          {loadError && <p>Run make site-export with completed pipeline trajectory results.</p>}
          {page === 'trajectories' && release?.assets.some((asset) => asset.kind === 'trajectory_preview') && (
            <figure className="trajectory-map-loading">
              <img
                src={`${import.meta.env.BASE_URL}data/pipeline-trajectory-preview.svg?release=${encodeURIComponent(release.releaseId)}`}
                alt="Preview of significant age-bin trajectories in the cross-cell-type map"
              />
              <figcaption role="status">Preparing the interactive map and its trajectory records…</figcaption>
            </figure>
          )}
        </section>
      )}
      {page === 'trajectories' && trajectoryMap && (
        <Explorer
          atlas={trajectoryMap}
          details={trajectoryDetails}
          onEnsureDetails={(point: TrajectoryPoint) => ensureTrajectoryDetails(point)}
          state={session.trajectory}
          preferredGene={session.selectedGene}
          onSelectGene={(gene, cellType) => dispatch({ type: 'selectGene', gene, cellType })}
          onStateChange={(state) => dispatch({ type: 'setTrajectory', state })}
          actions={pairActions}
          highlightGenes={session.highlightGenes}
          onHighlightGenesChange={updateHighlightGenes}
          pathwaySets={pathwaySets}
          activePathway={activePathway}
          onSelectPathway={selectPathway}
          geneMetadata={geneMetadata}
        />
      )}
      {page === 'genes' && atlas && (
        <ProfilesPage
          atlas={atlas}
          state={session.geneProfile}
          preferredGene={session.selectedGene}
          onSelectGene={(gene, cellType) => dispatch({ type: 'selectGene', gene, cellType })}
          onSelectCellType={(cellType) => dispatch({ type: 'selectCellType', cellType })}
          onStateChange={(state) => dispatch({ type: 'setProfile', tab: 'gene', state })}
          tab="gene"
          actions={pairActions}
          geneMetadata={geneMetadata}
        />
      )}
      {page === 'cellTypes' && atlas && (
        <ProfilesPage
          atlas={atlas}
          state={session.cellTypeProfile}
          preferredGene={session.selectedGene}
          onSelectGene={(gene, cellType) => dispatch({ type: 'selectGene', gene, cellType })}
          onSelectCellType={(cellType) => dispatch({ type: 'selectCellType', cellType })}
          onStateChange={(state) => dispatch({ type: 'setProfile', tab: 'cellType', state })}
          tab="cellType"
          actions={pairActions}
          highlightGenes={session.highlightGenes}
          onHighlightGenesChange={updateHighlightGenes}
          pathwaySets={pathwaySets}
          activePathway={activePathway}
          onSelectPathway={selectPathway}
          geneMetadata={geneMetadata}
        />
      )}
      {page === 'pathways' && (
        <PathwaySelector
          genes={session.highlightGenes}
          onGenesChange={updateHighlightGenes}
          pathwaySets={pathwaySets}
          activePathway={activePathway}
          onSelectPathway={selectPathway}
          onOpenMatrix={() => navigate('geneTypeMatrix')}
        />
      )}
      {page === 'geneTypeMatrix' && atlas && (
        <GeneTypeMatrix
          atlas={atlas}
          genes={session.highlightGenes}
          onGenesChange={updateHighlightGenes}
          onOpenGene={openGene}
          onOpenCellType={openCellType}
          pathwaySets={pathwaySets}
          activePathway={activePathway}
          onSelectPathway={selectPathway}
        />
      )}
      {page === 'composition' && (
        <PipelineComposition
          preferredCellType={session.selectedCellType}
          onSelectCellType={(cellType) => dispatch({ type: 'selectCellType', cellType })}
          denominator={session.compositionDenominator}
          onSelectDenominator={(denominator) =>
            dispatch({ type: 'selectCompositionDenominator', denominator })
          }
        />
      )}
      {page === 'singleCell' && <SingleCellPage preferredCellType={session.selectedCellType}
        onSelectCellType={(cellType) => dispatch({ type: 'selectCellType', cellType })} />}
      {(page === 'datasets' || page === 'modelling' || page === 'download') && (
        <InformationalPage page={page} />
      )}
      <footer>
        <span>Ageing Online Atlas</span>
        <span>Derived summary statistics · No clinical interpretation</span>
      </footer>
    </main>
  )
}
