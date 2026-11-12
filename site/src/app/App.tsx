import { useEffect, useReducer } from 'react'
import type { PathwaySet } from '../data/pathways'
import { CrossTypePage } from '../features/crossType/CrossTypePage'
import { DifferentialPage } from '../features/differential/DifferentialPage'
import { Explorer } from '../features/explorer/Explorer'
import { PipelineComposition } from '../features/composition/PipelineComposition'
import { InformationalPage } from '../features/content/InformationalPage'
import { OverviewPage } from '../features/content/OverviewPage'
import { PathwaysPage } from '../features/pathways/PathwaysPage'
import { ProfilesPage } from '../features/profiles/ProfilesPage'
import type { TrajectoryPoint, TrajectoryRecord } from '../data/atlas'
import { navigationItems } from './navigation'
import { NavigationMenu } from './NavigationMenu'
import { usePageNavigation } from './usePageNavigation'
import { useSiteData } from './useSiteData'
import { releaseLabel } from './releaseLabel'
import { atlasSessionReducer, initialAtlasSession } from './session'
import { SingleCellPage } from '../features/singleCell/SingleCellPage'
import { TranscriptomicClockPage } from '../features/content/TranscriptomicClockPage'
import { loadDifferential, preloadCombinedAgeModel } from '../data/differential'
import { initialCell } from '../features/differential/differentialView'

export function App() {
  const [session, dispatch] = useReducer(atlasSessionReducer, initialAtlasSession)
  useEffect(() => {
    let active = true
    loadDifferential()
      .then((data) => {
        if (!active) return
        const initial = initialCell(data.cellTypes, null, null)
        dispatch({ type: 'initializeCellType', cellType: initial.label })
        void preloadCombinedAgeModel(initial.label).catch(() => {})
      })
      .catch(() => {})
    return () => {
      active = false
    }
  }, [])
  const { page, navigate } = usePageNavigation()
  const {
    release,
    releaseLoading,
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
  const syntheticRelease =
    release?.analysis?.syntheticTestData === true || release?.inventory?.some((item) => item.syntheticTestData === true)
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
  const openProfile = (record: TrajectoryRecord, selectPair = true) => {
    dispatch({ type: 'selectGene', gene: record.gene, cellType: record.cellType })
    dispatch({ type: 'openProfile', tab: 'gene', record, selectPair })
    navigate('genes')
  }
  const openGene = (gene: string) => {
    dispatch({ type: 'selectGene', gene, cellType: session.selectedCellType ?? undefined })
    dispatch({
      type: 'setProfile',
      tab: 'gene',
      state: { ...session.geneProfile, gene, cellType: session.selectedCellType, selectedId: null },
    })
    navigate('genes')
  }
  const openGenePair = (gene: string, cellType: string) => {
    dispatch({ type: 'selectGene', gene, cellType })
    dispatch({ type: 'setProfile', tab: 'gene', state: { ...session.geneProfile, gene, cellType, selectedId: null } })
    navigate('genes')
  }
  const openCellType = (cellType: string) => {
    dispatch({ type: 'selectCellType', cellType })
    navigate('differential')
  }
  const openTrajectory = (record: TrajectoryRecord) => {
    dispatch({ type: 'selectGene', gene: record.gene, cellType: record.cellType })
    dispatch({ type: 'openTrajectory', record })
    navigate('trajectories')
  }
  const pairActions = {
    onOpenGene: (record: TrajectoryRecord) => openProfile(record),
    onOpenCellType: (record: TrajectoryRecord) => {
      dispatch({ type: 'selectGene', gene: record.gene, cellType: record.cellType })
      navigate('differential')
    },
    onOpenTrajectory: openTrajectory,
  }

  return (
    <main className={`atlas-main ${page === 'trajectories' ? 'is-trajectory-page' : ''}`}>
      <header className="site-header">
        <button className="brand" onClick={() => navigate('overview')}>
          <span className="brand-mark">B</span>
          <span>
            Blood Ageing
            <br />
            Explorer
          </span>
        </button>
        <nav className="site-nav" aria-label="Primary navigation">
          {navigationItems.map((item) => {
            if ('pages' in item) {
              return <NavigationMenu key={item.id} label={item.label} pages={item.pages} page={page} navigate={navigate} />
            }
            return (
              <button
                key={item.id}
                onClick={() => navigate(item.id)}
                className={page === item.id || (item.id === 'pathways' && page === 'geneTypeMatrix') ? 'active' : ''}
              >
                {item.label}
              </button>
            )
          })}
        </nav>
      </header>
      {annotationError &&
        annotationError !== 'unavailable' &&
        ['genes', 'trajectories'].includes(page) && (
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
      {page === 'overview' && (
        <OverviewPage onNavigate={navigate} summary={release?.context?.singleCellSummary} />
      )}
      {page === 'crossType' && <CrossTypePage preferredGene={session.selectedGene}
        onSelectGene={(gene) => dispatch({ type: 'selectGene', gene })} onOpenGene={openGene}
        onOpenCellType={(gene, cellType) => { dispatch({ type: 'selectGene', gene, cellType }); navigate('differential') }}
        onOpenTrajectories={() => navigate('trajectories')} />}
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
          onOpenGene={openGenePair}
          onOpenTrajectory={(gene, cellType) => {
            dispatch({ type: 'focusTrajectoryPair', gene, cellType })
            navigate('trajectories')
          }}
        />
      )}
      {!trajectoryMap && page === 'trajectories' && (
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
          {loadError && <p>Trajectory data could not be loaded for this release.</p>}
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
          preferredCellType={session.selectedCellType}
          onSelectGene={(gene, cellType) => dispatch({ type: 'selectGene', gene, cellType })}
          onStateChange={(state) => dispatch({ type: 'setTrajectory', state })}
          actions={pairActions}
          highlightGenes={session.highlightGenes}
          onHighlightGenesChange={updateHighlightGenes}
          pathwaySets={pathwaySets}
          activePathway={activePathway}
          onSelectPathway={selectPathway}
        />
      )}
      {page === 'genes' && (
        <ProfilesPage
          state={session.geneProfile}
          preferredGene={session.selectedGene}
          onSelectGene={(gene, cellType) => dispatch({ type: 'selectGene', gene, cellType })}
          onOpenCellType={openCellType}
          onOpenTrajectory={(gene, cellType) => {
            dispatch({ type: 'focusTrajectoryPair', gene, cellType })
            navigate('trajectories')
          }}
          onStateChange={(state) => dispatch({ type: 'setProfile', tab: 'gene', state })}
          geneMetadata={geneMetadata}
        />
      )}
      {['pathways', 'geneTypeMatrix'].includes(page) && (
        <PathwaysPage
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
      {(page === 'datasets' || page === 'modelling') && (
        <InformationalPage page={page} onNavigate={navigate} />
      )}
      {page === 'clock' && (
        <TranscriptomicClockPage
          preferredCellType={session.selectedCellType}
          onSelectCellType={(cellType) => dispatch({ type: 'selectCellType', cellType })}
        />
      )}
      {page === 'about' && <InformationalPage page="about" />}
      <footer>
        <div className="footer-branding">
          <span>Blood Ageing Explorer</span>
          <span>Derived summary statistics · No clinical interpretation</span>
          {syntheticRelease && <span>Synthetic test data · Not biological evidence</span>}
        </div>
        <span className="release" title={release ? `Export ID: ${release.releaseId}` : undefined}>
          {releaseLoading ? 'Loading release…' : releaseLabel(release)}
        </span>
      </footer>
    </main>
  )
}
