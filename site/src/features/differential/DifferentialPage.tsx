import { cellTypeColor } from '../../theme/cellTypeColors'
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
import { ChoiceMenu } from '../../components/ChoiceMenu'
import { VolcanoPlot } from '../../components/VolcanoPlot'
import { VolcanoExplorer } from '../../components/VolcanoExplorer'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import { choosePoint } from '../../data/pointSelection'
import { parseGeneSet } from '../../data/geneSet'
import type { PathwaySet } from '../../data/pathways'
import { visualColors } from '../../theme/visualEncoding'
import { defaultModelId, groupModels, modelLabel } from './modelSelection'
import { filterGenes, hasPlottableEstimate, initialCell, visibleGenes } from './differentialView'
import { studySupportColor } from './studySupportScale'
import { GeneAgePreview } from './GeneAgePreview'
import './differential.css'

const fmt = (value: number | null | undefined) => (value == null ? '—' : value.toPrecision(3))
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
  preferredCellType: string | null
  onSelectCellType: (cellType: string) => void
  preferredModelId: string | null
  onSelectModel: (modelId: string) => void
  pathwaySets: PathwaySet[]
  onSelectPathway: (pathway: PathwaySet) => void
  onHighlightGenesChange: (value: string) => void
}) {
  const entrySelection = useRef({ preferredGene, preferredCellType, preferredModelId, onSelectModel })
  const [release, setRelease] = useState<Release | null>(null)
  const [failure, setFailure] = useState<AssetFailureKind | null>(null)
  const [cellId, setCellId] = useState('')
  const [search, setSearch] = useState('')
  const [loadedRecords, setLoadedRecords] = useState<Record<string, Gene[]>>({})
  const [modelFailure, setModelFailure] = useState<AssetFailureKind | null>(null)
  const [lookupFailure, setLookupFailure] = useState(false)
  const [referenceExpression, setReferenceExpression] = useState<ReferenceExpression | null>(null)
  useEffect(() => {
    let active = true
    loadDifferential()
      .then((data) => {
        if (!active) return
        setRelease(data)
        const { preferredGene: gene, preferredCellType: type } = entrySelection.current
        const first = initialCell(data.cellTypes, gene, type)
        setCellId(first.id)
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
  const supportMaximum = points.reduce(
    (maximum, record) => Math.max(maximum, record.studySupport?.tested ?? 0),
    Math.max(1, Object.keys(model?.studies ?? {}).length),
  )
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
      {release.cellTypes.some((item) => item.models.some((entry) => !entry.recordsAsset)) && (
        <p role="status">
          This older export contains all models in one download. Re-export the site data to load models only
          when selected.
        </p>
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
            <p className="eyebrow">Model selection</p>{' '}
            <ChoiceMenu
              label="Cell type"
              value={cellId}
              groups={groupCellTypes([...new Set(release.cellTypes.map((item) => item.label))]).map(
                (group) => ({
                  label: group.label,
                  options: group.cellTypes.flatMap((label) =>
                    release.cellTypes
                      .filter((item) => item.label === label)
                      .map((item) => ({ value: item.id, label: item.label, color: cellTypeColor(item.label) })),
                  ),
                }),
              )}
              onChange={selectCell}
            />
            <ChoiceMenu
              label="Model and contrast"
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
            {modelId && !model && (
              <p role="status">
                The selected model is unavailable for {cell.label}. Select another model or choose a cell type
                with this result.
              </p>
            )}
            <label>
              Find a gene
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
        description={
          model ? (
            <>
              <p>
                Effect: {model.covariate.replaceAll('_', ' ')} {model.contrast && `· ${model.contrast}`}.
                Showing adjusted p &lt; 0.05 and explicitly requested genes. Enter a complete gene symbol to
                retrieve other results.
              </p>
              <p>Point color counts studies with adjusted p &lt; 0.05 for the same contrast, regardless of effect direction. Each gene's tooltip gives significant and tested study counts.</p>
              {geneSet.size > 0 && (
                <p>
                  Outlined: {highlightedCount.toLocaleString()} genes in{' '}
                  {activePathway?.label || 'the shared gene set'} with plottable estimates.
                </p>
              )}
            </>
          ) : (
            <p>Select an available model to view its results.</p>
          )
        }
        legend={
          <>
            {Array.from({ length: supportMaximum + 1 }, (_, count) => <span key={count}>
              <i style={{ backgroundColor: studySupportColor(count, supportMaximum) }} />
              {count}{count === supportMaximum ? ' studies' : ''}
            </span>)}
            <span><i style={{ backgroundColor: missingSupportColor }} /> Unavailable</span>
            <span>
              <i className="volcano-legend-pathway" /> Gene set
            </span>
            <span>
              <i className="volcano-legend-selected" /> Selected
            </span>
          </>
        }
        plot={
          model ? (
            <VolcanoPlot
              records={points}
              selected={selected}
              onSelect={(record) => onSelectGene(record.gene, cell.label)}
              id={(record) => record.gene}
              effect={(record) => record.effect}
              adjustedP={(record) => record.padj}
              highlighted={(record) => geneSet.has(record.gene.toUpperCase())}
              title={(record) =>
                `${record.gene}: effect ${fmt(record.effect)}, adjusted p ${fmt(record.padj)}; ${supportLabel(record.studySupport)}`
              }
              xLabel="Effect estimate (log₂ fold change)"
              threshold={0.05}
              ariaLabel="Volcano plot of effect versus negative log ten adjusted p-value. Search for a gene to select it with a keyboard."
              appearance={(record, isSelected) => ({
                fill: supportColor(record.studySupport, supportMaximum),
                stroke: isSelected
                  ? visualColors.selected
                  : geneSet.has(record.gene.toUpperCase())
                    ? visualColors.geneSetOutline
                    : 'none',
                strokeWidth: isSelected ? 2.5 : geneSet.has(record.gene.toUpperCase()) ? 1.5 : 0,
                radius: isSelected ? 7.5 : geneSet.has(record.gene.toUpperCase()) ? 6 : 4,
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
                  <h2>{selected.gene}</h2>
                  <GeneAgePreview
                    gene={selected.gene}
                    cellType={cell.label}
                    asset={cell.trajectoryAsset || null}
                  />
                  <p className="eyebrow">Selected model estimate</p>
                  <dl className="de-statistics">
                    <dt>Log₂ fold change</dt>
                    <dd>{fmt(selected.effect)}</dd>
                    <dt>Adjusted p-value</dt>
                    <dd>{fmt(selected.padj)}</dd>
                    <dt>Study support</dt>
                    <dd>{supportLabel(selected.studySupport)}</dd>
                    <dt>Standard error</dt>
                    <dd>{fmt(selected.se)}</dd>
                    <dt>Reference-bin expression ({referenceExpression?.referenceBin ?? '20–30'})</dt>
                    <dd>
                      {referenceExpression?.genes[selected.gene]
                        ? `${new Intl.NumberFormat('en', { maximumSignificantDigits: 4 }).format(referenceExpression.genes[selected.gene].cpm)} CPM · ${referenceExpression.genes[selected.gene].samples} samples, ${referenceExpression.genes[selected.gene].studies} studies`
                        : 'Unavailable'}
                    </dd>
                  </dl>
                  <p>
                    {referenceExpression?.method ||
                      'Study-balanced mean of per-study mean log1p(CPM), back-transformed with expm1.'}{' '}
                    CPM is normalized expression, not an absolute molecule count.
                  </p>
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
                    .join(', ') ||
                    model.study ||
                    'Not recorded'}
                </dd>
                <dt>Contrast</dt>
                <dd>{model.contrast || model.covariate}</dd>
              </dl>
              {model.downloadAsset && model.model !== 'per_study' && (
                <p className="de-model-download">
                  <a
                    href={`${import.meta.env.BASE_URL}data/${model.downloadAsset}?release=${encodeURIComponent(release.releaseId || '')}`}
                    download={`${cell.id}-${model.model}-${model.covariate}.csv`}
                  >
                    Download complete combined-model results (CSV)
                  </a>
                  <span>Includes all tested genes, including nonsignificant results.</span>
                </p>
              )}
            </div>
          ) : (
            <p>Select a model to see its details.</p>
          )
        }
      />
    </section>
  )
}
