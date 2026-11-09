import { useEffect, useMemo, useState } from 'react'
import { GeneSetInput } from '../../components/GeneSetInput'
import { PathwayAttribution } from '../../components/PathwayAttribution'
import type { Atlas } from '../../data/atlas'
import { loadCombinedAgeEffect } from '../../data/differential'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import { parseGeneSet } from '../../data/geneSet'
import type { PathwaySet } from '../../data/pathways'
import { formatAxisValue, symmetricAxis } from '../../theme/axisScale'
import { visualColors, visualRgb } from '../../theme/visualEncoding'
import { formatNumber } from '../../theme/format'
import './pathways.css'

export function PathwaysPage({
  atlas,
  loadError,
  genes,
  onGenesChange,
  onOpenGene,
  onOpenCellType,
  pathwaySets,
  activePathway,
  onSelectPathway,
}: {
  atlas: Atlas | null
  loadError: 'unavailable' | 'invalid' | 'mismatch' | null
  genes: string
  onGenesChange: (value: string) => void
  onOpenGene: (gene: string) => void
  onOpenCellType: (cellType: string) => void
  pathwaySets: PathwaySet[]
  activePathway: PathwaySet | null
  onSelectPathway: (pathway: PathwaySet) => void
}) {
  const [hideGenesWithoutResults, setHideGenesWithoutResults] = useState(false)
  const [combinedAgeEffects, setCombinedAgeEffects] = useState<Record<string, number | null>>({})
  const [effectsLoading, setEffectsLoading] = useState(false)
  const cellTypes = useMemo(
    () => atlas ? groupCellTypes(atlas.cellTypes).flatMap((group) => group.cellTypes) : [],
    [atlas?.cellTypes],
  )
  const selectedGenes = useMemo(() => [...parseGeneSet(genes)], [genes])
  useEffect(() => {
    let active = true
    if (!selectedGenes.length || !cellTypes.length) {
      setCombinedAgeEffects({})
      setEffectsLoading(false)
      return () => { active = false }
    }
    setEffectsLoading(true)
    void Promise.all(
      selectedGenes.flatMap((gene) => cellTypes.map(async (cellType) => {
        const key = `${gene}\u0000${cellType}`
        try {
          const yearlyEffect = await loadCombinedAgeEffect(cellType, gene)
          return [key, yearlyEffect === null ? null : yearlyEffect * 10] as const
        } catch {
          return [key, null] as const
        }
      })),
    ).then((entries) => {
      if (!active) return
      setCombinedAgeEffects(Object.fromEntries(entries))
      setEffectsLoading(false)
    })
    return () => { active = false }
  }, [cellTypes, selectedGenes])
  const matrix = useMemo(() => {
    return selectedGenes.map((gene) =>
      cellTypes.map((cellType) => combinedAgeEffects[`${gene}\u0000${cellType}`]),
    )
  }, [cellTypes, combinedAgeEffects, selectedGenes])
  const visibleGeneIndices = hideGenesWithoutResults
    ? selectedGenes.map((_, index) => index).filter((index) => matrix[index].some((value) => value != null))
    : selectedGenes.map((_, index) => index)
  const values = visibleGeneIndices.flatMap((index) =>
    matrix[index].flatMap((value) => (value == null ? [] : [value])),
  )
  const scale = symmetricAxis(values)
  const scalePosition = (value: number) => ((value + scale.extent) / (2 * scale.extent)) * 100
  const colorForValue = (value: number) => {
    const intensity = Math.min(1, Math.abs(value) / scale.extent)
    const [red, green, blue] = value >= 0 ? visualRgb.positiveChange : visualRgb.negativeChange
    return `rgb(${Math.round(255 + (red - 255) * intensity)} ${Math.round(255 + (green - 255) * intensity)} ${Math.round(255 + (blue - 255) * intensity)})`
  }
  const formatResultValue = formatNumber

  return (
    <section className="pathways-page">
      <header className="pathways-heading">
        <p className="eyebrow">Pathways</p>
        <h1>Pathway genes across cell types</h1>
        <p>Choose a curated gene set or enter symbols to compare their combined age effects across cell types.</p>
      </header>
      <div className="pathways-workspace">
        <div className="pathways-input-row">
          <div className="pathways-controls atlas-card">
            <GeneSetInput
              label="Gene set"
              value={genes}
              onChange={onGenesChange}
              helper="Enter symbols separated by commas or whitespace. Matching is exact and case-insensitive."
              pathwaySets={pathwaySets}
              activePathway={activePathway}
              onSelectPathway={onSelectPathway}
              pathwayDisplay="inline"
              onClear={() => onGenesChange('')}
            />
            {pathwaySets.length > 0 && <PathwayAttribution activePathway={activePathway} />}
            <p className="pathway-count">
              {selectedGenes.length.toLocaleString()} {selectedGenes.length === 1 ? 'gene' : 'genes'} selected
            </p>
          </div>
          <div
            className="heatmap-scale"
            aria-label={
              selectedGenes.length > 0
                ? `Symmetric combined age-model effect scale from ${formatAxisValue(-scale.extent, scale.extent)} to ${formatAxisValue(scale.extent, scale.extent)} log2 fold change per decade`
                : undefined
            }
          >
            {selectedGenes.length > 0 && (
              <>
                <p className="heatmap-scale-title">Combined age model β · log₂ fold change per decade</p>
                {values.length > 0 ? (
                  <div className="heatmap-scale-ramp">
                    <div
                      className="heatmap-scale-bar"
                      style={{
                        background: `linear-gradient(to right, ${colorForValue(-scale.extent)} 0%, #fffefa 50%, ${colorForValue(scale.extent)} 100%)`,
                      }}
                    />
                    {scale.ticks.map((tick) => {
                      const position = scalePosition(tick)
                      return (
                        <span
                          className="heatmap-scale-label"
                          key={tick}
                          style={{
                            left: `${position}%`,
                            transform:
                              position <= 0
                                ? undefined
                                : position >= 100
                                  ? 'translateX(-100%)'
                                  : 'translateX(-50%)',
                          }}
                        >
                          {formatAxisValue(tick, scale.extent)}
                        </span>
                      )
                    })}
                  </div>
                ) : effectsLoading ? (
                  <p className="pathway-no-matches" role="status">Loading combined age effects…</p>
                ) : (
                  <p className="pathway-no-matches">
                    No combined age-model estimates are retained for these genes.
                  </p>
                )}
                <label className="pathway-filter">
                  <input
                    type="checkbox"
                    checked={hideGenesWithoutResults}
                    onChange={(event) => setHideGenesWithoutResults(event.target.checked)}
                  />{' '}
                  Hide genes without combined age estimates
                </label>
              </>
            )}
          </div>
        </div>
        <div className="pathways-matrix-heading">
          <div>
            <h2>Combined age effect</h2>
            <p>Each value is the combined linear-age model estimate in log₂ fold change per decade. Select a gene or cell type label to open its detailed view.</p>
          </div>
        </div>
        {!atlas ? (
          <div className="empty-pathways" role={loadError ? 'alert' : 'status'}>
            {loadError === 'invalid'
              ? 'The age-effect data is invalid.'
              : loadError === 'mismatch'
                ? 'The age-effect data belongs to a different release.'
                : loadError
                  ? 'Age-effect data is unavailable for this release.'
                  : 'Loading age-effect data…'}
          </div>
        ) : selectedGenes.length === 0 ? (
          <div className="empty-pathways">
            Choose a curated gene set or enter gene symbols above to display the matrix.
          </div>
        ) : visibleGeneIndices.length === 0 ? (
          <div className="empty-pathways">
            No selected genes have combined age-model estimates in this dataset.
          </div>
        ) : effectsLoading ? (
          <div className="empty-pathways" role="status">Loading combined age effects…</div>
        ) : (
          <div className="heatmap-scroll atlas-card atlas-card-soft">
            <div
              className="heatmap"
              style={{
                gridTemplateColumns: `110px repeat(${cellTypes.length}, 32px)`,
                gridTemplateRows: `105px repeat(${visibleGeneIndices.length}, 32px)`,
              }}
            >
              <div className="heatmap-corner">Gene / cell type</div>
              {cellTypes.map((cellType) => (
                <button
                  className="heatmap-column-label"
                  type="button"
                  key={cellType}
                  onClick={() => onOpenCellType(cellType)}
                >
                  <span>{cellType}</span>
                </button>
              ))}
              {visibleGeneIndices.flatMap((geneIndex) => {
                const gene = selectedGenes[geneIndex]
                return [
                  <button
                    className="heatmap-row-label"
                    type="button"
                    key={`${gene}-label`}
                    onClick={() => onOpenGene(gene)}
                  >
                    {gene}
                  </button>,
                  ...matrix[geneIndex].map((record, cellIndex) => (
                    <div
                      className="heatmap-cell"
                      key={`${gene}-${cellTypes[cellIndex]}`}
                      style={{ backgroundColor: record != null ? colorForValue(record) : visualColors.missingResult }}
                      title={record != null ? `${record > 0 ? '+' : ''}${formatResultValue(record)} log₂FC per decade` : 'No retained combined age-model estimate'}
                    >
                      <span className="sr-only">
                        {record != null
                          ? `${gene}, ${cellTypes[cellIndex]}: ${record > 0 ? '+' : ''}${formatResultValue(record)} log2 fold change per decade`
                          : `${gene}, ${cellTypes[cellIndex]}: no retained combined age-model estimate`}
                      </span>
                    </div>
                  )),
                ]
              })}
            </div>
          </div>
        )}
      </div>
    </section>
  )
}
