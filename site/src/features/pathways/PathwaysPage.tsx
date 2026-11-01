import { useMemo, useState } from 'react'
import { GeneSetInput } from '../../components/GeneSetInput'
import { PathwayAttribution } from '../../components/PathwayAttribution'
import type { Atlas, TrajectoryRecord } from '../../data/atlas'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import { parseGeneSet } from '../../data/geneSet'
import type { PathwaySet } from '../../data/pathways'
import { formatAxisValue, symmetricAxis } from '../../theme/axisScale'
import { visualColors, visualRgb } from '../../theme/visualEncoding'
import './pathways.css'

function recordWithLargestAbsoluteChange(records: TrajectoryRecord[]): TrajectoryRecord | undefined {
  return records.reduce<TrajectoryRecord | undefined>(
    (largest, record) =>
      !largest || Math.abs(record.valueAtMaxAbsChange) > Math.abs(largest.valueAtMaxAbsChange)
        ? record
        : largest,
    undefined,
  )
}

export function GeneTypeMatrix({
  atlas,
  genes,
  onGenesChange,
  onOpenGene,
  onOpenCellType,
  pathwaySets,
  activePathway,
  onSelectPathway,
}: {
  atlas: Atlas
  genes: string
  onGenesChange: (value: string) => void
  onOpenGene: (gene: string) => void
  onOpenCellType: (cellType: string) => void
  pathwaySets: PathwaySet[]
  activePathway: PathwaySet | null
  onSelectPathway: (pathway: PathwaySet) => void
}) {
  const [hideGenesWithoutResults, setHideGenesWithoutResults] = useState(false)
  const cellTypes = useMemo(
    () => groupCellTypes(atlas.cellTypes).flatMap((group) => group.cellTypes),
    [atlas.cellTypes],
  )
  const selectedGenes = useMemo(() => [...parseGeneSet(genes)], [genes])
  const matrix = useMemo(() => {
    const selectedGeneSet = new Set(selectedGenes)
    const recordsByGeneAndCellType = new Map<string, TrajectoryRecord[]>()
    atlas.records.forEach((record) => {
      const normalizedGene = record.gene.toUpperCase()
      if (!selectedGeneSet.has(normalizedGene)) return
      const key = `${normalizedGene}\u0000${record.cellType}`
      recordsByGeneAndCellType.set(key, [...(recordsByGeneAndCellType.get(key) ?? []), record])
    })
    return selectedGenes.map((gene) =>
      cellTypes.map((cellType) =>
        recordWithLargestAbsoluteChange(recordsByGeneAndCellType.get(`${gene}\u0000${cellType}`) ?? []),
      ),
    )
  }, [atlas.records, cellTypes, selectedGenes])
  const visibleGeneIndices = hideGenesWithoutResults
    ? selectedGenes.map((_, index) => index).filter((index) => matrix[index].some(Boolean))
    : selectedGenes.map((_, index) => index)
  const values = visibleGeneIndices.flatMap((index) =>
    matrix[index].flatMap((record) => (record ? [record.valueAtMaxAbsChange] : [])),
  )
  const scale = symmetricAxis(values)
  const scalePosition = (value: number) => ((value + scale.extent) / (2 * scale.extent)) * 100
  // The scroll frame has a 1px border on both sides; include it in the workspace
  // width so a grid that exactly fits does not create a 2px horizontal overflow.
  const heatmapWidth = 112 + cellTypes.length * 32
  const colorForValue = (value: number) => {
    const intensity = Math.min(1, Math.abs(value) / scale.extent)
    const [red, green, blue] = value >= 0 ? visualRgb.positiveChange : visualRgb.negativeChange
    return `rgb(${Math.round(255 + (red - 255) * intensity)} ${Math.round(255 + (green - 255) * intensity)} ${Math.round(255 + (blue - 255) * intensity)})`
  }
  const color = (record: TrajectoryRecord | undefined) =>
    record ? colorForValue(record.valueAtMaxAbsChange) : visualColors.missingResult
  const formatResultValue = (value: number) => value.toFixed(2)

  return (
    <section className="pathways-page">
      <h1 className="visually-hidden">Gene × cell type</h1>
      <div className="pathways-workspace" style={{ width: `min(100%, max(760px, ${heatmapWidth}px))` }}>
        <div className="pathways-input-row">
          <div className="pathways-controls">
            <GeneSetInput
              label="Highlight Genes"
              value={genes}
              onChange={onGenesChange}
              helper="Separate symbols with commas or whitespace. This selection is synchronized to the other data tabs. Changes here will be reflected there as well."
              pathwaySets={pathwaySets}
              activePathway={activePathway}
              onSelectPathway={onSelectPathway}
              pathwayDisplay="inline"
            />
            {pathwaySets.length > 0 && <PathwayAttribution activePathway={activePathway} />}
          </div>
          <div
            className="heatmap-scale"
            aria-label={
              selectedGenes.length > 0
                ? `Symmetric signed maximum absolute age change scale from ${formatAxisValue(-scale.extent, scale.extent)} to ${formatAxisValue(scale.extent, scale.extent)}`
                : undefined
            }
          >
            {selectedGenes.length > 0 && (
              <>
                <p className="heatmap-scale-title">Maximum change from 20–30-year-old baseline</p>
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
                ) : (
                  <p className="pathway-no-matches">
                    No selected genes have retained trajectory results in this dataset.
                  </p>
                )}
                <label className="pathway-filter">
                  <input
                    type="checkbox"
                    checked={hideGenesWithoutResults}
                    onChange={(event) => setHideGenesWithoutResults(event.target.checked)}
                  />{' '}
                  Hide genes without retained trajectory results
                </label>
              </>
            )}
          </div>
        </div>
        {selectedGenes.length === 0 ? (
          <div className="empty-pathways">
            Paste one or more retained gene symbols to display the heatmap.
          </div>
        ) : visibleGeneIndices.length === 0 ? (
          <div className="empty-pathways">
            No selected genes have retained trajectory results in this dataset.
          </div>
        ) : (
          <div className="heatmap-scroll">
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
                      style={{ backgroundColor: color(record) }}
                      title={record ? formatResultValue(record.valueAtMaxAbsChange) : 'No retained result'}
                    >
                      <span className="sr-only">
                        {record
                          ? `${gene}, ${cellTypes[cellIndex]}: ${formatResultValue(record.valueAtMaxAbsChange)}`
                          : `${gene}, ${cellTypes[cellIndex]}: no retained result`}
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
