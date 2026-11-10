import { assetFailureKind, type AssetFailureKind } from '../../data/release'
import {
  loadDifferential,
  loadModelRecords,
  loadReferenceExpression,
  type ReferenceExpression,
} from '../../data/differential'
import { geneBucket } from '../../data/geneBucket.mjs'
import type { DifferentialRelease as Release, Gene } from '../../data/contracts.mjs'
import { useEffect, useMemo, useRef, useState } from 'react'
import { CellTypeSelector } from '../../components/CellTypeSelector'
import { ChoiceMenu } from '../../components/ChoiceMenu'
import { VolcanoPlot } from '../../components/VolcanoPlot'
import { VolcanoExplorer } from '../../components/VolcanoExplorer'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import { choosePoint } from '../../data/pointSelection'
import { parseGeneSet } from '../../data/geneSet'
import type { PathwaySet } from '../../data/pathways'
import { visualColors } from '../../theme/visualEncoding'
import { formatCpm, formatNumber } from '../../theme/format'
import { defaultModelId, groupModels, modelLabel } from './modelSelection'
import { filterGenes, hasPlottableEstimate, initialCell, visibleGenes } from './differentialView'
import { studySupportColor, viridisColor } from './studySupportScale'
import { DifferentialDownload } from './DifferentialDownload'
import { GeneAgePreview } from './GeneAgePreview'
import './differential.css'

const missingSupportColor = '#b9bdbc'
const supportColor = (support: Gene['studySupport'], maximum: number) =>
  support ? studySupportColor(support.significant, maximum) : missingSupportColor
const supportLabel = (support: Gene['studySupport']) =>
  support ? `${support.significant} of ${support.tested} studies significant` : 'Study support unavailable'

export function DifferentialPage({
  loadingPreviewPath,
  releaseId,
  highlightGenes,
  activePathway,
  preferredGene,
  onSelectGene,
  onOpenGene,
  onOpenTrajectory,
  preferredCellType,
  onSelectCellType,
  preferredModelId,
  onSelectModel,
  pathwaySets,
  onSelectPathway,
  onHighlightGenesChange,
}: {
  loadingPreviewPath?: string
  releaseId?: string
  highlightGenes: string
  activePathway: PathwaySet | null
  preferredGene: string | null
  onSelectGene: (gene: string | null, cellType?: string) => void
  onOpenGene: (gene: string, cellType: string) => void
  onOpenTrajectory: (gene: string, cellType: string) => void
  preferredCellType: string | null
  onSelectCellType: (cellType: string) => void
  preferredModelId: string | null
  onSelectModel: (modelId: string) => void
  pathwaySets: PathwaySet[]
  onSelectPathway: (pathway: PathwaySet) => void
  onHighlightGenesChange: (value: string) => void
}) {
  const entrySelection = useRef({ preferredGene, preferredCellType, preferredModelId, onSelectModel, onSelectCellType })
  const [release, setRelease] = useState<Release | null>(null)
  const [failure, setFailure] = useState<AssetFailureKind | null>(null)
  const [cellId, setCellId] = useState('')
  const [search, setSearch] = useState('')
  const [loadedRecords, setLoadedRecords] = useState<Record<string, Gene[]>>({})
  const [modelFailure, setModelFailure] = useState<AssetFailureKind | null>(null)
  const [lookupFailure, setLookupFailure] = useState(false)
  const [referenceExpression, setReferenceExpression] = useState<ReferenceExpression | null>(null)
  const [colorPreference, setColorPreference] = useState<{ key: string; enabled: boolean } | null>(null)
  useEffect(() => {
    let active = true
    loadDifferential()
      .then((data) => {
        if (!active) return
        setRelease(data)
        const { preferredGene: gene, preferredCellType: type } = entrySelection.current
        const first = initialCell(data.cellTypes, gene, type)
        setCellId(first.id)
        entrySelection.current.onSelectCellType(first.label)
        const defaultId = defaultModelId(first.models)
        if (entrySelection.current.preferredModelId === null && defaultId)
          entrySelection.current.onSelectModel(defaultId)
      })
      .catch((error) => {
        if (active) setFailure(assetFailureKind(error))
      })
    return () => {
      active = false
    }
  }, [])
  const cell = release?.cellTypes.find((item) => item.id === cellId)
  const modelId = preferredModelId ?? (cell ? defaultModelId(cell.models) : '')
  const modelMetadata = cell?.models.find((item) => item.id === modelId)
  const ageEffectScale = modelMetadata?.covariate === 'age' && modelMetadata.contrast === 'per year' ? 10 : 1
  const recordsAsset = modelMetadata?.recordsAsset
  useEffect(() => {
    if (!cell?.referenceExpressionAsset) {
      setReferenceExpression(null)
      return
    }
    let active = true
    setReferenceExpression(null)
    loadReferenceExpression(cell.referenceExpressionAsset, cell.label)
      .then((value) => {
        if (active) setReferenceExpression(value)
      })
      .catch(() => {
        if (active) setReferenceExpression(null)
      })
    return () => {
      active = false
    }
  }, [cell?.referenceExpressionAsset, cell?.label])
  useEffect(() => {
    if (!recordsAsset || !cell || !modelMetadata || loadedRecords[recordsAsset]) return
    let active = true
    setModelFailure(null)
    loadModelRecords(recordsAsset, cell.id, modelMetadata.id)
      .then((records) => {
        if (active) setLoadedRecords((previous) => ({ ...previous, [recordsAsset]: records }))
      })
      .catch((error) => {
        if (active) setModelFailure(assetFailureKind(error))
      })
    return () => {
      active = false
    }
  }, [recordsAsset, cell, modelMetadata, loadedRecords])
  const requestedGenes = useMemo(
    () => new Set([search.trim().toUpperCase(), preferredGene?.toUpperCase() || ''].filter(Boolean)),
    [search, preferredGene],
  )
  useEffect(() => {
    if (!cell || !modelMetadata?.recordBuckets) return
    const paths = [
      ...new Set(
        [...requestedGenes]
          .map((gene) => modelMetadata.recordBuckets?.[geneBucket(gene)])
          .filter((path): path is string => !!path && !loadedRecords[path]),
      ),
    ]
    if (!paths.length) return
    let active = true
    const timer = window.setTimeout(() => {
      setLookupFailure(false)
      Promise.all(
        paths.map(async (path) => [path, await loadModelRecords(path, cell.id, modelMetadata.id)] as const),
      )
        .then((results) => {
          if (active) setLoadedRecords((previous) => ({ ...previous, ...Object.fromEntries(results) }))
        })
        .catch(() => {
          if (active) setLookupFailure(true)
        })
    }, 250)
    return () => {
      active = false
      window.clearTimeout(timer)
    }
  }, [cell, modelMetadata, loadedRecords, requestedGenes])
  useEffect(() => {
    if (!cell || !modelMetadata?.recordBuckets) return
    const paths = [
      ...new Set(
        [...parseGeneSet(highlightGenes)]
          .map((gene) => modelMetadata.recordBuckets?.[geneBucket(gene)])
          .filter((path): path is string => !!path),
      ),
    ]
    let active = true
    let cursor = 0
    const worker = async () => {
      while (active && cursor < paths.length) {
        const path = paths[cursor++]
        try {
          await loadModelRecords(path, cell.id, modelMetadata.id)
        } catch {
          /* Explicit lookups report failures. */
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
  }, [cell, modelMetadata, highlightGenes])
  const model = modelMetadata && {
    ...modelMetadata,
    records: visibleGenes(
      recordsAsset
        ? [
            ...(loadedRecords[recordsAsset] || []),
            ...[...requestedGenes].flatMap((gene) => {
              const path = modelMetadata.recordBuckets?.[geneBucket(gene)]
              return path
                ? (loadedRecords[path] || []).filter((record) => record.gene.toUpperCase() === gene)
                : []
            }),
          ]
        : modelMetadata.records,
      requestedGenes,
    ),
  }
  const requestedModel = release?.cellTypes.flatMap((item) => item.models).find((item) => item.id === modelId)
  const filtered = useMemo(() => filterGenes(model, search), [model, search])
  const geneSet = useMemo(() => parseGeneSet(highlightGenes), [highlightGenes])
  if (!release || !cell)
    return (
      <section className="content-page">
        <p className="eyebrow">Differential expression</p>
        <h1>
          {failure === 'invalid'
            ? 'Differential expression data is invalid.'
            : failure === 'mismatch'
              ? 'Differential expression data belongs to a different release.'
              : failure
                ? 'No differential expression export yet.'
                : 'Loading differential expression…'}
        </h1>
        {failure && (
          <p>
            Run <code>make site-export</code> with completed pipeline results to populate this view.
          </p>
        )}
        {loadingPreviewPath && (
          <figure className="de-preview">
            <img
              src={`${import.meta.env.BASE_URL}data/${loadingPreviewPath}?release=${encodeURIComponent(releaseId || '')}`}
              alt="Static preview of the combined age differential-expression results"
            />
            <figcaption role="status">Loading the interactive differential-expression results…</figcaption>
          </figure>
        )}
      </section>
    )
  const points = (model?.records || []).filter(hasPlottableEstimate)
  const selected = choosePoint(
    points,
    (record) => record.padj,
    (record) => record.gene === preferredGene,
  )
  const highlightedCount = points.filter((record) => geneSet.has(record.gene.toUpperCase())).length
  const hasStudySupport = points.some((record) => record.studySupport != null)
  const supportMaximum = Math.max(
    1,
    ...points.flatMap((record) => record.studySupport ? [record.studySupport.tested] : []),
  )
  const referenceExpressionValues = Object.values(referenceExpression?.genes ?? {})
    .filter((gene) => Number.isFinite(gene.cpm) && gene.cpm >= 0)
    .map((gene) => Math.log2(gene.cpm + 1))
  const expressionMinimum = referenceExpressionValues.length ? Math.min(...referenceExpressionValues) : 0
  const expressionMaximum = referenceExpressionValues.length ? Math.max(...referenceExpressionValues) : 0
  const colorPreferenceKey = `${cellId}:${modelId}`
  const colorByReferenceExpression = colorPreference?.key === colorPreferenceKey
    ? colorPreference.enabled
    : !hasStudySupport
  const referenceExpressionColor = (gene: string) => {
    const cpm = referenceExpression?.genes[gene]?.cpm
    if (cpm == null || !Number.isFinite(cpm) || cpm < 0) return missingSupportColor
    const value = Math.log2(cpm + 1)
    const normalized = expressionMaximum > expressionMinimum
      ? (value - expressionMinimum) / (expressionMaximum - expressionMinimum)
      : 0.5
    return viridisColor(normalized)
  }
  const pointColorLabel = (record: Gene) =>
    colorByReferenceExpression
      ? referenceExpression?.genes[record.gene]
        ? `log₂(reference-bin CPM + 1): ${formatNumber(Math.log2(referenceExpression.genes[record.gene].cpm + 1))}`
        : 'Reference expression unavailable'
      : hasStudySupport
        ? supportLabel(record.studySupport)
        : 'Study support unavailable'
  const combinedDownload =
    model && (model.model === 'merged' || model.model === 'combined') && model.downloadAsset
      ? model.downloadAsset
      : null
  const selectCell = (id: string) => {
    const next = release.cellTypes.find((item) => item.id === id)!
    setCellId(id)
    onSelectCellType(next.label)
  }
  return (
    <section className="de-page">
      <h1 className="visually-hidden">Differential expression</h1>
      {preferredCellType && !release.cellTypes.some((item) => item.label === preferredCellType) && (
        <p role="status">
          Differential expression results are unavailable for {preferredCellType}. Showing {cell.label}.
        </p>
      )}
      {(cell.warning || cell.analysisMode === 'test_only') && (
        <div className="de-warning" role="status">
          {cell.warning || 'Test-only results are not biological evidence.'}
        </div>
      )}
      {recordsAsset && !loadedRecords[recordsAsset] && (
        <p role="status">
          {modelFailure ? 'Could not load the selected model results.' : 'Loading selected model results…'}
        </p>
      )}
      {modelMetadata?.model === 'merged' && modelMetadata.covariate === 'age' &&
        recordsAsset && !loadedRecords[recordsAsset] && cell.previewAsset && (
          <figure className="de-preview">
            <img
              src={`${import.meta.env.BASE_URL}data/${cell.previewAsset}?release=${encodeURIComponent(release?.releaseId || releaseId || '')}`}
              alt={`Static preview of the combined age differential-expression results for ${cell.label}`}
            />
            <figcaption role="status">Loading the interactive results for {cell.label}…</figcaption>
          </figure>
        )}
      {lookupFailure && <p role="status">Could not load requested gene results.</p>}
      <VolcanoExplorer
        controls={
          <>
            <p className="eyebrow">Explore data</p>
            <CellTypeSelector
              label="Cell type"
              value={cellId}
              groups={groupCellTypes([...new Set(release.cellTypes.map((item) => item.label))]).map(
                (group) => ({
                  label: group.label,
                  options: group.cellTypes.flatMap((label) =>
                    release.cellTypes
                      .filter((item) => item.label === label)
                      .map((item) => ({ value: item.id, label: item.label })),
                  ),
                }),
              )}
              onChange={selectCell}
            />
            <ChoiceMenu
              label="Analysis"
              value={modelId}
              displayValue={
                modelId && !model
                  ? `${requestedModel ? modelLabel(requestedModel) : modelId} (unavailable)`
                  : undefined
              }
              groups={groupModels(cell.models).map((group) => ({
                label: group.label,
                options: group.models.map((item) => ({ value: item.id, label: modelLabel(item) })),
              }))}
              onChange={onSelectModel}
            />
            <div className="de-combined-download">
              {combinedDownload ? (
                <DifferentialDownload
                  key={`${release.releaseId}:${combinedDownload}`}
                  path={combinedDownload}
                  releaseId={release.releaseId || ''}
                  filename={`${cell.id}-${model?.covariate}-combined.csv`}
                />
              ) : (
                <button className="de-download-button" type="button" disabled>
                  Download DEG results
                </button>
              )}
              <p>
                {combinedDownload
                  ? `Combined ${model?.covariate} results for ${cell.label}. The CSV starts with # metadata lines; skip them to read the table in software.`
                  : model && (model.model === 'merged' || model.model === 'combined')
                    ? 'A download is not available for this result in the current release.'
                    : 'Downloads are available for combined models only.'}
              </p>
            </div>
            {modelId && !model && (
              <p role="status">
                The selected model is unavailable for {cell.label}. Select another model or choose a cell type
                with this result.
              </p>
            )}
            <label>
              Search genes
              <input
                type="search"
                value={search}
                onChange={(event) => setSearch(event.target.value)}
                placeholder="Gene symbol"
              />
            </label>
            {model && search.trim() && (
              <div className="de-gene-matches">
                <ChoiceMenu
                  label={`Matching genes (${filtered.length.toLocaleString()})`}
                  value={selected && filtered.includes(selected) ? selected.gene : ''}
                  displayValue={selected && filtered.includes(selected) ? undefined : 'Select a gene'}
                  groups={[
                    {
                      options: [
                        { value: '', label: 'Clear selection' },
                        ...filtered.slice(0, 100).map((gene) => ({ value: gene.gene, label: gene.gene })),
                      ],
                    },
                  ]}
                  onChange={(value) => onSelectGene(value || null, cell.label)}
                />
                {filtered.length > 100 && <span>Showing the first 100 matches; refine your search.</span>}
              </div>
            )}
          </>
        }
        geneSet={{
          value: highlightGenes,
          onChange: onHighlightGenesChange,
          pathwaySets,
          activePathway,
          onSelectPathway,
        }}
        title="Volcano plot"
        showTitle={false}
        descriptionPlacement="below"
        legendPlacement="below"
        description={
          model ? (
            <>
              <p>
                Effect: {model.covariate.replaceAll('_', ' ')} {model.contrast && `· ${model.contrast}`}.
                Showing adjusted p &lt; 0.05 and explicitly requested genes. Enter a complete gene symbol to
                retrieve other results.
              </p>
              <p>
                {colorByReferenceExpression
                  ? 'Point color shows log₂(reference-bin CPM + 1). Genes without a reference expression value are gray.'
                  : hasStudySupport
                    ? 'Point color counts studies with adjusted p < 0.05 for the same contrast, regardless of direction. Each gene’s tooltip gives significant and tested study counts.'
                    : 'Study support is unavailable for this model, so points are gray.'}
              </p>
              {geneSet.size > 0 && (
                <p>
                  Outlined: {highlightedCount.toLocaleString()} genes in {activePathway?.label || 'the shared gene set'}
                  {' '}with plottable estimates.
                </p>
              )}
            </>
          ) : (
            <p>Select an available model to view its results.</p>
          )
        }
        legend={
          <>
            {colorByReferenceExpression ? (
              <>
                <span className="de-support-legend" aria-label="Log base two reference-bin CPM plus one">
                  <span className="de-support-label">log₂(reference-bin CPM + 1)</span>
                  <span className="de-support-scale">
                    <span className="de-support-gradient" aria-hidden="true" />
                    <span className="de-support-endpoints">
                      <span>{formatNumber(expressionMinimum)}</span>
                      <span aria-hidden="true">–</span>
                      <span>{formatNumber(expressionMaximum)}</span>
                    </span>
                  </span>
                </span>
                {points.some((record) => referenceExpression?.genes[record.gene] == null) && (
                  <span className="de-support-unavailable"><i />Not available</span>
                )}
              </>
            ) : hasStudySupport ? (
              <>
                <span className="de-support-legend" aria-label="Number of studies with a significant result">
                  <span className="de-support-label">Significant studies</span>
                  <span className="de-support-scale">
                    <span className="de-support-gradient" aria-hidden="true" />
                    <span className="de-support-endpoints"><span>1</span><span aria-hidden="true">–</span><span>{supportMaximum}</span></span>
                  </span>
                </span>
                <span className="de-support-zero"><i />0</span>
                {points.some((record) => record.studySupport == null) && (
                  <span className="de-support-unavailable"><i />Not available</span>
                )}
              </>
            ) : <span className="de-support-unavailable"><i />Study support unavailable</span>}
            <span>
              <i className="volcano-legend-pathway" /> Gene set
            </span>
            <span>
              <i className="volcano-legend-selected" /> Selected
            </span>
          </>
        }
        plot={
          model && recordsAsset && !loadedRecords[recordsAsset] ? (
            <p className="de-plot-loading" role="status">
              {modelFailure ? 'The selected model results could not be loaded.' : 'Loading selected model results…'}
            </p>
          ) : model ? (
            <VolcanoPlot
              records={points}
              selected={selected}
              onSelect={(record) => onSelectGene(record.gene, cell.label)}
              id={(record) => record.gene}
              effect={(record) => record.effect * ageEffectScale}
              adjustedP={(record) => record.padj}
              highlighted={(record) => geneSet.has(record.gene.toUpperCase())}
              aspectRatio="square"
              className="de-volcano-plot"
              title={(record) =>
                `${record.gene}: effect ${formatNumber(record.effect * ageEffectScale)}${ageEffectScale === 10 ? ' log₂ fold change per decade' : ''}, adjusted p ${formatNumber(record.padj)}; ${pointColorLabel(record)}`
              }
              xLabel={ageEffectScale === 10 ? 'Combined age effect estimates (log₂ fold change per decade)' : 'Effect estimate (log₂ fold change)'}
              threshold={0.05}
              ariaLabel="Volcano plot of effect versus negative log ten adjusted p-value. Search for a gene to select it with a keyboard."
              appearance={(record, isSelected) => ({
                fill: colorByReferenceExpression
                  ? referenceExpressionColor(record.gene)
                  : hasStudySupport
                    ? supportColor(record.studySupport, supportMaximum)
                    : missingSupportColor,
                stroke: isSelected
                  ? visualColors.selected
                  : geneSet.has(record.gene.toUpperCase())
                    ? visualColors.geneSetOutline
                    : 'none',
                strokeWidth: isSelected
                  ? 2.5
                  : geneSet.has(record.gene.toUpperCase())
                    ? 2.5
                    : 0,
                radius: isSelected ? 9.5 : geneSet.has(record.gene.toUpperCase()) ? 7 : 5.5,
                opacity:
                  !isSelected && search && !record.gene.toUpperCase().includes(search.trim().toUpperCase())
                    ? 0.12
                    : 0.8,
              })}
            />
          ) : (
            <p>No model selected.</p>
          )
        }
        detail={
          model ? (
            <div className="de-context">
              <p className="eyebrow">Selected gene</p>
              {preferredGene && !points.some((record) => record.gene === preferredGene) && (
                <p>
                  The selected gene {preferredGene} has no plottable estimate in this model. Showing another
                  result.
                </p>
              )}
              {selected ? (
                <>
                  <h2>
                    <button className="de-gene-link" onClick={() => onOpenGene(selected.gene, cell.label)}>
                      {selected.gene}
                    </button>
                  </h2>
                  {model.covariate === 'age' && cell.trajectoryAsset && (
                    <GeneAgePreview
                      gene={selected.gene}
                      cellType={cell.label}
                      asset={cell.trajectoryAsset}
                    />
                  )}
                  <p className="eyebrow">Selected model estimate</p>
                  <dl className="de-statistics">
                    <dt>{ageEffectScale === 10 ? 'Log₂ fold change per decade' : 'Log₂ fold change'}</dt>
                    <dd>{formatNumber(selected.effect * ageEffectScale)}</dd>
                    <dt>Adjusted p-value</dt>
                    <dd>{formatNumber(selected.padj)}</dd>
                    <dt>Study support</dt>
                    <dd>{supportLabel(selected.studySupport)}</dd>
                    <dt>{ageEffectScale === 10 ? 'Standard error per decade' : 'Standard error'}</dt>
                    <dd>{formatNumber(selected.se === null ? null : selected.se * ageEffectScale)}</dd>
                    <dt>Reference-bin expression (CPM)</dt>
                    <dd>
                      {referenceExpression?.genes[selected.gene]
                        ? formatCpm(referenceExpression.genes[selected.gene].cpm)
                        : 'Unavailable'}
                      <label className="de-expression-toggle" title="Color points by reference expression">
                        <input
                          type="checkbox"
                          aria-label="Color points by reference expression"
                          checked={colorByReferenceExpression}
                          disabled={!referenceExpression}
                          onChange={(event) =>
                            setColorPreference({ key: colorPreferenceKey, enabled: event.target.checked })
                          }
                        />
                        <span>Color</span>
                      </label>
                    </dd>
                  </dl>
                  <button
                    className="text-button"
                    onClick={() => onOpenTrajectory(selected.gene, cell.label)}
                  >
                    View in age-bin trajectories →
                  </button>
                </>
              ) : (
                <p>Select a point or search for a gene to inspect its estimate.</p>
              )}
              <p className="eyebrow">Model context</p>
              <h2>{modelLabel(model)}</h2>
              <dl>
                <dt>Design</dt>
                <dd>{model.design}</dd>
                <dt>Samples</dt>
                <dd>{model.samples}</dd>
                <dt>Studies</dt>
                <dd>
                  {Object.entries(model.studies)
                    .map(([name, count]) => `${name} (n=${count})`)
                    .join(', ') || model.study || 'Not recorded'}
                </dd>
                <dt>Contrast</dt>
                <dd>{model.contrast || model.covariate}</dd>
              </dl>
            </div>
          ) : (
            <p>Select a model to see its details.</p>
          )
        }
      />
    </section>
  )
}
