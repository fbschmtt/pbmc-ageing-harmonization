import { lazy, Suspense, useEffect, useMemo } from 'react'
import { GeneTypePairCard, type GeneTypePairActions } from '../../components/GeneTypePairCard'
import { GeneSetInput } from '../../components/GeneSetInput'
import type { TrajectoryMap, TrajectoryPoint, TrajectoryRecord } from '../../data/atlas'
import { choosePoint } from '../../data/pointSelection'
import { parseGeneSet } from '../../data/geneSet'
import type { PathwaySet } from '../../data/pathways'
import './explorer.css'
import { compileFocusQuery, hasRegexSpecialCharacters, recordMatchesFocus } from './explorerFilter'
import type { UmapViewState } from './umapView'

const Umap = lazy(() => import('./Umap').then((module) => ({ default: module.Umap })))

export type ExplorerState = {
  query: string
  selectedId: string | null
  umapViewState: UmapViewState | null
  umapResetVersion: number
}
type ExplorerProps = {
  atlas: TrajectoryMap
  details: TrajectoryRecord[]
  onEnsureDetails: (point: TrajectoryPoint) => void
  state: ExplorerState
  onStateChange: (state: ExplorerState) => void
  preferredGene: string | null
  preferredCellType: string | null
  onSelectGene: (gene: string, cellType?: string) => void
  actions: GeneTypePairActions
  highlightGenes: string
  onHighlightGenesChange: (value: string) => void
  pathwaySets: PathwaySet[]
  activePathway: PathwaySet | null
  onSelectPathway: (pathway: PathwaySet) => void
}
export function Explorer({
  atlas,
  details,
  onEnsureDetails,
  state,
  onStateChange,
  preferredGene,
  preferredCellType,
  onSelectGene,
  actions,
  highlightGenes,
  onHighlightGenesChange,
  pathwaySets,
  activePathway,
  onSelectPathway,
}: ExplorerProps) {
  const updateState = (update: Partial<ExplorerState>) => onStateChange({ ...state, ...update })
  const regexResult = useMemo(() => compileFocusQuery(state.query), [state.query])
  const filtered = useMemo(
    () => atlas.records.filter((record) => recordMatchesFocus(record, state.query, regexResult.expression)),
    [atlas.records, regexResult.expression, state.query],
  )
  const selectedInView =
    filtered.find((record) => record.id === state.selectedId) ??
    choosePoint(
      filtered,
      (record) => record.adjustedPValue,
      (record) => record.gene === preferredGene && record.cellType === preferredCellType,
    )
  const highlightedGenes = useMemo(() => parseGeneSet(highlightGenes), [highlightGenes])
  const selectedDetail = details.find((record) => record.id === selectedInView?.id)
  useEffect(() => {
    if (selectedInView) onEnsureDetails(selectedInView)
  }, [selectedInView?.id, onEnsureDetails])
  const highlightedMatchCount = useMemo(
    () => filtered.filter((record) => highlightedGenes.has(record.gene.toUpperCase())).length,
    [filtered, highlightedGenes],
  )
  return (
    <section className="explorer">
      <h1 className="explorer-heading">Age-bin trajectories</h1>
      <div className="explorer-layout">
        <aside className="filters">
          <label>
            Focus by gene or cell type
            <input
              value={state.query}
              onChange={(event) => updateState({ query: event.target.value, selectedId: null })}
              placeholder="e.g. PDK4, OPC, or ^Oligo"
            />
          </label>
          <button
            className="filter-clear"
            type="button"
            disabled={!state.query}
            onClick={() => updateState({ query: '', selectedId: null })}
          >
            Clear focus
          </button>
          <button
            className="map-reset"
            type="button"
            onClick={() => updateState({ umapViewState: null, umapResetVersion: state.umapResetVersion + 1 })}
          >
            Center UMAP
          </button>
          {regexResult.error && (
            <p className="search-error" role="alert">
              {regexResult.error}
            </p>
          )}
          {hasRegexSpecialCharacters(state.query) && (
            <p className="search-warning">Special characters are interpreted as a regular expression.</p>
          )}
          <GeneSetInput
            className="trajectory-highlight-genes"
            label="Highlight gene set"
            value={highlightGenes}
            onChange={onHighlightGenesChange}
            onClear={() => onHighlightGenesChange('')}
            helper="Dark, outlined points belong to this set. The set is synchronized with the Pathways tab."
            pathwaySets={pathwaySets}
            activePathway={activePathway}
            onSelectPathway={onSelectPathway}
          />
          <p className="focus-status" aria-live="polite">
            {state.query
              ? `${filtered.length.toLocaleString()} focus matches`
              : `${atlas.records.length.toLocaleString()} points`}
            {highlightedGenes.size > 0 &&
              ` · ${highlightedMatchCount.toLocaleString()} ${state.query ? 'also ' : ''}in gene set`}
          </p>
        </aside>
        <div className="map-panel">
          <Suspense fallback={<div className="empty-map">Loading the similarity map…</div>}>
            <Umap
              records={atlas.records}
              matchingRecords={filtered}
              selected={selectedInView}
              focusSelection={Boolean(state.selectedId)}
              highlightedGenes={highlightedGenes}
              initialViewState={state.umapViewState}
              resetVersion={state.umapResetVersion}
              onViewStateChange={(umapViewState) => updateState({ umapViewState })}
              onSelect={(record) => {
                onSelectGene(record.gene, record.cellType)
                updateState({ selectedId: record.id })
              }}
            />
          </Suspense>
          <div className="legend">
            <span>
              <i className="brain-dot" /> Brain
            </span>
            <span>
              <i className="blood-dot" /> Peripheral blood
            </span>
            <span>
              <i className="gene-set-ring" /> Gene set
            </span>
            <span>
              <i className="selected-ring" /> Selected
            </span>
          </div>
        </div>
        <aside className="detail-panel">
          {selectedDetail ? (
            <GeneTypePairCard
              record={selectedDetail}
              actions={{ ...actions, onOpenTrajectory: undefined }}
            />
          ) : selectedInView ? (
            <p role="status">Loading details for {selectedInView.gene}…</p>
          ) : (
            <p>
              {filtered.length
                ? 'Select a point for additional information.'
                : 'No results match the current focus.'}
            </p>
          )}
        </aside>
      </div>
    </section>
  )
}
