import { useEffect, useMemo, useState } from 'react'
import { VolcanoExplorer } from '../../components/VolcanoExplorer'
import { GeneFunctionSummary } from '../../components/GeneFunctionSummary'
import { loadDifferential, loadModelRecords, loadCombinedAgeResult } from '../../data/differential'
import type { DifferentialRelease, Gene } from '../../data/contracts.mjs'
import type { GeneMetadataIndex } from '../../data/geneMetadata'
import { choosePoint } from '../../data/pointSelection'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import { formatNumber } from '../../theme/format'
import { ProfileScatter, type CombinedAgePoint } from './ProfileScatter'
import { GeneAgePreview } from '../differential/GeneAgePreview'
import './profiles.css'
import '../../components/geneTypePairCard.css'

export type ProfileTab = 'gene' | 'cellType'
export type ProfileState = { gene: string | null; cellType: string | null; selectedId: string | null }

export function ProfilesPage({ state, onStateChange, preferredGene, onSelectGene, onOpenCellType,
  onOpenTrajectory, geneMetadata = null }: {
  state: ProfileState
  onStateChange: (state: ProfileState) => void
  preferredGene: string | null
  onSelectGene: (gene: string, cellType?: string) => void
  onOpenCellType: (cellType: string) => void
  onOpenTrajectory: (gene: string, cellType: string) => void
  geneMetadata?: GeneMetadataIndex | null
}) {
  const [catalogue, setCatalogue] = useState<DifferentialRelease | null>(null)
  const [defaultGene, setDefaultGene] = useState('')
  const [options, setOptions] = useState<string[]>([])
  const [loaded, setLoaded] = useState<{ gene: string; records: CombinedAgePoint[] } | null>(null)
  const [failed, setFailed] = useState(false)
  const gene = (state.gene ?? preferredGene ?? defaultGene).trim().toUpperCase()
  useEffect(() => {
    let active = true
    loadDifferential().then((value) => { if (active) setCatalogue(value) })
      .catch(() => { if (active) setFailed(true) })
    return () => { active = false }
  }, [])
  const cells = useMemo(() => catalogue?.cellTypes.filter((cell) =>
    cell.models.some((model) => model.model === 'merged' && model.covariate === 'age')) ?? [], [catalogue])
  useEffect(() => {
    if (!cells.length || gene) return
    let active = true
    Promise.all(cells.map((cell) => {
      const model = cell.models.find((model) => model.model === 'merged' && model.covariate === 'age')!
      return model.recordsAsset ? loadModelRecords(model.recordsAsset, cell.id, model.id) : model.records
    })).then((rows) => {
      if (!active) return
      const records = rows.flat()
      setOptions([...new Set(records.map((row) => row.gene))].sort())
      setDefaultGene(choosePoint(records, (row: Gene) => row.padj)?.gene ?? records[0]?.gene ?? '')
    }).catch(() => { if (active) setFailed(true) })
    return () => { active = false }
  }, [cells, gene])
  useEffect(() => {
    if (!catalogue || !gene) return
    let active = true
    setFailed(false)
    const timer = window.setTimeout(() => {
      Promise.all(cells.map(async (cell) => {
        const row = await loadCombinedAgeResult(cell.label, gene)
        if (row?.effect == null || row.padj == null) return null
        return { ...row, effect: row.effect, padj: row.padj,
          id: `${row.gene} | ${cell.label}`, cellType: cell.label, cellId: cell.id,
          trajectoryAsset: cell.trajectoryAsset } satisfies CombinedAgePoint
      })).then((rows) => {
        if (active) setLoaded({ gene, records: rows.filter((row): row is NonNullable<typeof row> => row !== null) })
      }).catch(() => { if (active) setFailed(true) })
    }, 250)
    return () => { active = false; window.clearTimeout(timer) }
  }, [catalogue, cells, gene])
  const records = loaded?.gene === gene ? loaded.records : []
  const selected = records.find((row) => row.id === state.selectedId) ??
    choosePoint(records, (row) => row.padj, (row) => row.cellType === state.cellType)
  const ordered = groupCellTypes(records.map((row) => row.cellType)).flatMap((group) => group.cellTypes)
  const sortedRecords = [...records].sort((a, b) => ordered.indexOf(a.cellType) - ordered.indexOf(b.cellType))
  return <section className="profiles-page">
    <h1 className="visually-hidden">Combined age effects by gene</h1>
    <VolcanoExplorer showTitle={false} descriptionPlacement="below" legendPlacement="below" title="Combined age effects by cell type"
      description={<p>Each point represents the combined linear age model for this gene in one cell type. Effects are log₂ fold changes per decade; adjusted p-values come from the same model.</p>}
      legend={<><span>Point color: AIFI L2 cell type (gray for other types)</span><span>Selected point is outlined</span></>}
      controls={<>
        <p className="eyebrow">Gene selection</p>
        <label>Gene symbol<input list="gene-options" value={state.gene ?? preferredGene ?? defaultGene}
          onChange={(event) => {
            const value = event.target.value
            if (value.trim()) onSelectGene(value.trim().toUpperCase())
            onStateChange({ ...state, gene: value, selectedId: null })
          }} /></label>
        <datalist id="gene-options">{options.map((item) => <option key={item} value={item} />)}</datalist>
        <GeneFunctionSummary metadata={geneMetadata?.genes[gene]} release={geneMetadata?.release} />
        <p>Compare this gene across cell types with completed combined age models.</p>
      </>}
      plot={failed ? <p role="alert">Combined age-model results could not be loaded.</p>
        : !catalogue || (gene && loaded?.gene !== gene) ? <p role="status">Loading combined age-model results…</p>
        : !records.length ? <p>No finite combined age-model estimates are available for this gene.</p>
        : <ProfileScatter records={sortedRecords} selected={selected}
          onSelect={(row) => {
            onSelectGene(row.gene, row.cellType)
            onStateChange({ ...state, gene: row.gene, cellType: row.cellType, selectedId: row.id })
          }} />}
      detail={selected ? <article className="pair-card">
        <p className="eyebrow">Gene–celltype pair</p><h2>{selected.gene}</h2>
        <p className="detail-subtitle"><button className="pair-card-cell-link" onClick={() => onOpenCellType(selected.cellType)}>{selected.cellType}</button> · pbmc</p>
        <dl><div><dt>Combined linear age effect per decade</dt><dd>{formatNumber(selected.effect * 10)} log₂ fold change</dd></div>
          <div><dt>Adjusted p-value (combined age model)</dt><dd>{formatNumber(selected.padj)}</dd></div>
          <div><dt>Standard error per decade</dt><dd>{formatNumber(selected.se === null ? null : selected.se * 10)}</dd></div></dl>
        <GeneAgePreview gene={selected.gene} cellType={selected.cellType} asset={selected.trajectoryAsset} />
        <div className="pair-card-actions">
          <button className="text-button" onClick={() => onOpenCellType(selected.cellType)}>View differential expression →</button>
          <button className="text-button" onClick={() => onOpenTrajectory(selected.gene, selected.cellType)}>Open age-bin trajectory map →</button>
        </div>
      </article> : <p>Select a point for additional information.</p>}
    />
  </section>
}
