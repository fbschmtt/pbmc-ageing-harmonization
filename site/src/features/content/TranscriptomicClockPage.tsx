import { useEffect, useMemo, useState } from 'react'
import { CellTypeSelector } from '../../components/CellTypeSelector'
import { ChoiceMenu } from '../../components/ChoiceMenu'
import { loadClockIndex, loadClockModel } from '../../data/clocks'
import type { ClockCellType, ClockMethod, ClockMethodSummary, ClockPrediction } from '../../data/contracts.mjs'
import { formatNumber } from '../../theme/format'
import './transcriptomicClocks.css'

const palette = ['#287d70', '#c5795d', '#7896a5', '#a18a4e', '#775f8d', '#598b62', '#b46469']
const format = formatNumber

function PredictionPlot({
  title,
  points,
}: {
  title: string
  points: ClockPrediction[]
}) {
  const studies = [...new Set(points.map((point) => point.study))].sort()
  const colors = new Map(studies.map((study, index) => [study, palette[index % palette.length]]))
  const min = 0, max = 100
  const range = max - min
  const width = 560, left = 58, right = 18, top = 18, bottom = 50
  const height = top + (width - left - right) + bottom
  const x = (value: number) => left + ((value - min) / range) * (width - left - right)
  const y = (value: number) => height - bottom - ((value - min) / range) * (height - top - bottom)
  const clippedX = (value: number) => x(Math.max(min, Math.min(max, value)))
  const clippedY = (value: number) => y(Math.max(min, Math.min(max, value)))
  const ticks = Array.from({ length: 5 }, (_, index) => min + (range * index) / 4)
  const outsideCount = points.filter((point) => point.age < min || point.age > max || point.predictedAge < min || point.predictedAge > max).length
  return (
    <section className="clock-plot" aria-label={title}>
      <h3>{title}</h3>
      <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={`${title}: observed age compared with predicted age`}>
        {ticks.map((tick) => (
          <g key={tick}>
            <line x1={left} x2={width - right} y1={y(tick)} y2={y(tick)} className="clock-grid" />
            <line x1={x(tick)} x2={x(tick)} y1={top} y2={height - bottom} className="clock-grid" />
            <text x={x(tick)} y={height - bottom + 20} textAnchor="middle" className="clock-tick">{Math.round(tick)}</text>
            <text x={left - 10} y={y(tick) + 4} textAnchor="end" className="clock-tick">{Math.round(tick)}</text>
          </g>
        ))}
        <line x1={x(min)} y1={y(min)} x2={x(max)} y2={y(max)} className="clock-reference" />
        <g>
          {points.map((point, index) => {
            const outsideX = point.age < min || point.age > max
            const outsideY = point.predictedAge < min || point.predictedAge > max
            const label = `${point.study}: observed ${format(point.age)} years, predicted ${format(point.predictedAge)} years${point.nCellTypes === undefined ? '' : `; mean of ${point.nCellTypes} cell types`}`
            if (!outsideX && !outsideY) return (
              <circle key={`${point.study}-${point.age}-${index}`} cx={x(point.age)} cy={y(point.predictedAge)} r="4" fill={colors.get(point.study)} className="clock-dot">
                <title>{label}</title>
              </circle>
            )
            const dx = outsideX ? Math.sign(point.age - (min + max) / 2) : 0
            const dy = outsideY ? -Math.sign(point.predictedAge - (min + max) / 2) : 0
            const angle = Math.atan2(dy, dx) * 180 / Math.PI
            return (
              <path key={`${point.study}-${point.age}-${index}`} d="M 6 0 L -4 -4 L -2 0 L -4 4 Z" transform={`translate(${clippedX(point.age)} ${clippedY(point.predictedAge)}) rotate(${angle})`} fill={colors.get(point.study)} className="clock-outlier-arrow">
                <title>{`${label}; outside the 0–100 year plot range`}</title>
              </path>
            )
          })}
        </g>
        <text x={(left + width - right) / 2} y={height - 8} textAnchor="middle" className="clock-axis">Observed age (years)</text>
        <text transform={`translate(15 ${(top + height - bottom) / 2}) rotate(-90)`} textAnchor="middle" className="clock-axis">Predicted age (years)</text>
      </svg>
      <p className="clock-outlier-note" aria-hidden={outsideCount === 0}>
        {outsideCount > 0 ? `${outsideCount} ${outsideCount === 1 ? 'point is' : 'points are'} beyond 0–100 years; hover an arrow for its value.` : '\u00a0'}
      </p>
      <div className="clock-study-key" aria-label="Studies">
        {studies.map((study) => <span key={study}><i style={{ background: colors.get(study) }} />{study}</span>)}
      </div>
    </section>
  )
}

function Performance({ metrics }: { metrics: { n: number; mae: number; rmse: number; r2: number | null } }) {
  return <p className="clock-metrics">MAE {format(metrics.mae)} years <span>·</span> RMSE {format(metrics.rmse)} years <span>·</span> R² {metrics.r2 === null ? '—' : format(metrics.r2)} <span>·</span> n = {metrics.n}</p>
}

function ClockCard({ method }: { method: ClockMethod }) {
  const isCrossType = method.id.startsWith('cross_type_mean:')
  return (
    <article className="clock-method-card atlas-card">
      <header>
        <p className="eyebrow">{isCrossType ? 'Cross-cell-type ensemble' : method.id === 'elastic_net' ? 'Baseline clock' : 'Clock model'}</p>
        <h2>{method.name}</h2>
        <p>{method.description}</p>
      </header>
      <div className="clock-plots">
        <div>
          <Performance metrics={method.training.metrics} />
          <PredictionPlot title="Training studies · leave-one-study-out" points={method.training.predictions} />
        </div>
        <div>
          <Performance metrics={method.heldOut.metrics} />
          <PredictionPlot title="Held-out studies" points={method.heldOut.predictions} />
        </div>
      </div>
    </article>
  )
}

export function TranscriptomicClockPage({
  preferredCellType,
  onSelectCellType,
}: {
  preferredCellType: string | null
  onSelectCellType: (cellType: string) => void
}) {
  const [cellTypes, setCellTypes] = useState<ClockCellType[]>([])
  const [selectedCellId, setSelectedCellId] = useState('')
  const [selectedMethodId, setSelectedMethodId] = useState('')
  const [loadedMethod, setLoadedMethod] = useState<ClockMethod | null>(null)
  const [status, setStatus] = useState<'loading' | 'ready' | 'unavailable' | 'invalid'>('loading')
  const [modelStatus, setModelStatus] = useState<'loading' | 'ready' | 'unavailable' | 'invalid'>('loading')

  useEffect(() => {
    let active = true
    loadClockIndex().then((index) => {
      if (!active) return
      setCellTypes(index.cellTypes)
      setStatus('ready')
    }).catch((error) => {
      if (!active) return
      setStatus(error?.kind === 'invalid' ? 'invalid' : 'unavailable')
    })
    return () => { active = false }
  }, [])

  useEffect(() => {
    if (!cellTypes.length) return
    const preferred = cellTypes.find((cell) => cell.scope !== 'cross_cell_type' && cell.label === preferredCellType)
    const firstCell = cellTypes.find((cell) => cell.scope !== 'cross_cell_type') ?? cellTypes[0]
    setSelectedCellId((preferred ?? firstCell).id)
  }, [cellTypes, preferredCellType])

  const selectedCell = useMemo(() => cellTypes.find((cell) => cell.id === selectedCellId), [cellTypes, selectedCellId])
  const selectedSummary: ClockMethodSummary | undefined = selectedCell?.methods.find((method) => method.id === selectedMethodId)

  useEffect(() => {
    if (!selectedCell?.methods.length) return
    if (!selectedCell.methods.some((method) => method.id === selectedMethodId))
      setSelectedMethodId(selectedCell.methods[selectedCell.methods.length - 1].id)
  }, [selectedCell, selectedMethodId])

  useEffect(() => {
    if (!selectedCell || !selectedSummary) return
    let active = true
    setLoadedMethod(null)
    setModelStatus('loading')
    loadClockModel(selectedSummary.asset, selectedCell.id, selectedSummary.id).then((model) => {
      if (!active) return
      setLoadedMethod({
        id: selectedSummary.id,
        name: selectedSummary.name,
        description: selectedSummary.description,
        training: { ...selectedSummary.training, predictions: model.training.predictions },
        heldOut: { ...selectedSummary.heldOut, predictions: model.heldOut.predictions },
      })
      setModelStatus('ready')
    }).catch((error) => {
      if (!active) return
      setModelStatus(error?.kind === 'invalid' ? 'invalid' : 'unavailable')
    })
    return () => { active = false }
  }, [selectedCell, selectedSummary])

  const selectCellType = (id: string) => {
    setSelectedCellId(id)
    const cell = cellTypes.find((item) => item.id === id)
    if (cell && cell.scope !== 'cross_cell_type') onSelectCellType(cell.label)
  }

  return (
    <section className="content-page clock-page">
      <p className="eyebrow">Transcriptomic clocks</p>
      <h1>Predicting age from gene expression</h1>
      <p className="clock-intro">Compare expression-based clocks by cell type and model. Training performance uses leave-one-study-out predictions; held-out performance comes from studies excluded from model fitting and tuning. “Across cell types” shows the donor-matched mean prediction when available.</p>
      {status === 'loading' && <p role="status">Loading clock results…</p>}
      {status === 'unavailable' && <p role="status">Clock results are not available in this site release.</p>}
      {status === 'invalid' && <p role="alert">Clock results could not be read because the release data are invalid.</p>}
      {status === 'ready' && cellTypes.length === 0 && <p>No cell types have complete clock results in this release.</p>}
      {cellTypes.length > 0 && selectedCell && (
        <>
          <div className="clock-selectors">
            <CellTypeSelector
              className="clock-type-select"
              value={selectedCell.id}
              groups={[{
                options: cellTypes.map((cell) => ({
                  value: cell.id,
                  label: cell.label,
                  ...(cell.scope === 'cross_cell_type' ? { color: '#287d70' } : {}),
                })),
              }]}
              onChange={selectCellType}
            />
            <div className="clock-model-select">
              <ChoiceMenu
                label="Clock model"
                value={selectedSummary?.id ?? ''}
                groups={[{ options: selectedCell.methods.map((method) => ({ value: method.id, label: method.name })) }]}
                onChange={setSelectedMethodId}
              />
            </div>
          </div>
          <p className="clock-explanation">Each dot is a donor pseudobulk. Axes use the same 0–100 year range for every model and cell type; arrows mark points beyond that range. The diagonal marks perfect agreement between observed and predicted age. Study colors make differences in cohort performance visible.</p>
          {modelStatus === 'loading' && <p role="status">Loading selected clock model…</p>}
          {modelStatus === 'unavailable' && <p role="status">The selected clock model is not available in this site release.</p>}
          {modelStatus === 'invalid' && <p role="alert">The selected clock model could not be read because the release data are invalid.</p>}
          {modelStatus === 'ready' && loadedMethod && <ClockCard method={loadedMethod} />}
        </>
      )}
    </section>
  )
}
