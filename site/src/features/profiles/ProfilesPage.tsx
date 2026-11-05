import { cellTypeColor } from '../../theme/cellTypeColors'
import { useMemo } from 'react'
import { ChoiceMenu } from '../../components/ChoiceMenu'
import { GeneTypePairCard, type GeneTypePairActions } from '../../components/GeneTypePairCard'
import { VolcanoExplorer } from '../../components/VolcanoExplorer'
import { GeneFunctionSummary } from '../../components/GeneFunctionSummary'
import type { Atlas } from '../../data/atlas'
import type { GeneMetadataIndex } from '../../data/geneMetadata'
import { choosePoint } from '../../data/pointSelection'
import { parseGeneSet } from '../../data/geneSet'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import type { PathwaySet } from '../../data/pathways'
import { ProfileScatter } from './ProfileScatter'
import './profiles.css'

export type ProfileTab = 'gene' | 'cellType'
export type ProfileState = { gene: string | null; cellType: string | null; selectedId: string | null }

export function ProfilesPage({
  atlas,
  state,
  onStateChange,
  preferredGene,
  onSelectGene,
  onSelectCellType,
  tab,
  actions,
  highlightGenes = '',
  onHighlightGenesChange,
  pathwaySets = [],
  activePathway = null,
  onSelectPathway,
  geneMetadata = null,
}: {
  atlas: Atlas
  state: ProfileState
  onStateChange: (state: ProfileState) => void
  preferredGene: string | null
  onSelectGene: (gene: string, cellType?: string) => void
  onSelectCellType: (cellType: string) => void
  tab: ProfileTab
  actions: GeneTypePairActions
  highlightGenes?: string
  onHighlightGenesChange?: (value: string) => void
  pathwaySets?: PathwaySet[]
  activePathway?: PathwaySet | null
  onSelectPathway?: (pathway: PathwaySet) => void
  geneMetadata?: GeneMetadataIndex | null
}) {
  const defaultGene = choosePoint(atlas.records, (record) => record.adjustedPValue)?.gene ?? atlas.genes[0]
  const gene =
    state.gene ?? (preferredGene && atlas.genes.includes(preferredGene) ? preferredGene : defaultGene)
  const cellType = state.cellType ?? atlas.cellTypes[0]
  const parsedHighlightGenes = useMemo(() => parseGeneSet(highlightGenes), [highlightGenes])
  const cellTypeGroups = useMemo(() => groupCellTypes(atlas.cellTypes), [atlas.cellTypes])
  const records = useMemo(
    () =>
      atlas.records.filter((record) =>
        tab === 'gene' ? record.gene === gene : record.cellType === cellType,
      ),
    [atlas.records, tab, gene, cellType],
  )
  const selected =
    records.find((record) => record.id === state.selectedId) ??
    choosePoint(
      records,
      (record) => record.adjustedPValue,
      (record) => record.gene === preferredGene,
    )
  const contextualActions =
    tab === 'gene' ? { ...actions, onOpenGene: undefined } : { ...actions, onOpenCellType: undefined }
  return (
    <section className="profiles-page">
      <h1 className="visually-hidden">
        {tab === 'cellType' ? 'Explore trajectories by cell type' : 'Explore trajectories by gene'}
      </h1>
      <VolcanoExplorer
        title={tab === 'cellType' ? 'Age-bin effects by gene' : 'Age-bin effects by cell type'}
        description={
          <p>
            Each point represents a gene–cell type combination. Effects show the largest signed age-bin
            change.
          </p>
        }
        legend={
          <>
            {tab === 'gene' ? (
              <span>Point color: AIFI L2 cell type (gray for other types)</span>
            ) : (
              <>
                <span><i className="volcano-legend-blood" /> Peripheral blood</span>
                <span><i className="volcano-legend-brain" /> Brain</span>
              </>
            )}
            {tab === 'cellType' && (
              <span>
                <i className="volcano-legend-pathway" /> Gene set
              </span>
            )}
            <span>
              <i className="volcano-legend-selected" /> Selected
            </span>
          </>
        }
        geneSet={
          tab === 'cellType'
            ? {
                value: highlightGenes,
                onChange: (value) => onHighlightGenesChange?.(value),
                pathwaySets,
                activePathway,
                onSelectPathway,
              }
            : undefined
        }
        controls={
          <>
            <p className="eyebrow">{tab === 'gene' ? 'Gene selection' : 'Celltype selection'}</p>
            {tab === 'gene' ? (
              <label>
                Gene symbol
                <input
                  list="gene-options"
                  value={gene}
                  onChange={(event) => {
                    const value = event.target.value
                    if (atlas.genes.includes(value)) onSelectGene(value)
                    onStateChange({ ...state, gene: value, selectedId: null })
                  }}
                />
                <datalist id="gene-options">
                  {atlas.genes.map((item) => (
                    <option key={item} value={item} />
                  ))}
                </datalist>
              </label>
            ) : (
              <ChoiceMenu
                label="Cell type"
                value={cellType}
                groups={cellTypeGroups.map((group) => ({
                  label: group.label,
                  options: group.cellTypes.map((item) => ({ value: item, label: item, color: cellTypeColor(item) })),
                }))}
                onChange={(value) => {
                  onSelectCellType(value)
                  onStateChange({ ...state, cellType: value, selectedId: null })
                }}
              />
            )}
            {tab === 'gene' && (
              <GeneFunctionSummary metadata={geneMetadata?.genes[gene]} release={geneMetadata?.release} />
            )}
            {tab === 'gene' && <p>Compare this gene across all retained cell type and tissue results.</p>}
          </>
        }
        plot={
          <ProfileScatter
            records={records}
            selected={selected}
            onSelect={(record) => {
              onSelectGene(record.gene, record.cellType)
              onStateChange({ ...state, gene: record.gene, cellType: record.cellType, selectedId: record.id })
            }}
            label={(record) => (tab === 'gene' ? record.cellType : record.gene)}
            colorByCellType={tab === 'gene'}
            highlightedGenes={tab === 'cellType' ? parsedHighlightGenes : undefined}
          />
        }
        detail={
          selected ? (
            <GeneTypePairCard
              record={selected}
              actions={contextualActions}
            />
          ) : (
            <div className="empty-scatter">Select a point for additional information.</div>
          )
        }
      />
    </section>
  )
}
