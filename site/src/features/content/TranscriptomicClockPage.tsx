import { useEffect, useMemo, useState } from 'react'
import { loadClockIndex } from '../../data/clocks'
import type { ClockCellType, ClockMethod, ClockPrediction } from '../../data/contracts.mjs'
import './transcriptomicClocks.css'

const palette = ['#287d70', '#c5795d', '#7896a5', '#a18a4e', '#775f8d', '#598b62', '#b46469']
const format = (value: number) => value.toFixed(1)

function PredictionPlot({ title, points }: { title: string; points: ClockPrediction[] }) {
  const studies = [...new Set(points.map((point) => point.study))].sort()
  const colors = new Map(studies.map((study, index) => [study, palette[index % palette.length]]))
  const allValues = points.flatMap((point) => [point.age, point.predictedAge])
  const min = Math.floor(Math.min(...allValues) / 10) * 10
  const max = Math.ceil(Math.max(...allValues) / 10) * 10
  const range = Math.max(max - min, 1)
  const width = 560, height = 370, left = 58, right = 18, top = 18, bottom = 50
  const x = (value: number) => left + ((value - min) / range) * (width - left - right)
  const y = (value: number) => height - bottom - ((value - min) / range) * (height - top - bottom)
  const ticks = Array.from({ length: 5 }, (_, index) => min + (range * index) / 4)
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
        {points.map((point, index) => (
          <circle key={`${point.study}-${point.age}-${index}`} cx={x(point.age)} cy={y(point.predictedAge)} r="4" fill={colors.get(point.study)} className="clock-dot">
            <title>{`${point.study}: observed ${format(point.age)} years, predicted ${format(point.predictedAge)} years${point.nCellTypes === undefined ? '' : `; mean of ${point.nCellTypes} cell types`}`}</title>
          </circle>
        ))}
        <text x={(left + width - right) / 2} y={height - 8} textAnchor="middle" className="clock-axis">Observed age (years)</text>
        <text transform={`translate(15 ${(top + height - bottom) / 2}) rotate(-90)`} textAnchor="middle" className="clock-axis">Predicted age (years)</text>
      </svg>
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
  return (
    <article className="clock-method-card atlas-card">
      <header>
        <p className="eyebrow">{method.id.startsWith('cross_type_mean:') ? 'Cross-cell-type ensemble' : method.id === 'elastic_net' ? 'Baseline clock' : 'Most advanced available'}</p>
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

export function TranscriptomicClockPage() {
  const [cellTypes, setCellTypes] = useState<ClockCellType[]>([])
  const [selectedId, setSelectedId] = useState('')
  const [status, setStatus] = useState<'loading' | 'ready' | 'unavailable' | 'invalid'>('loading')
  useEffect(() => {
    let active = true
    loadClockIndex().then((index) => {
      if (!active) return
      setCellTypes(index.cellTypes)
      setSelectedId(index.cellTypes[0]?.id ?? '')
      setStatus('ready')
    }).catch((error) => {
      if (!active) return
      setStatus(error?.kind === 'invalid' ? 'invalid' : 'unavailable')
    })
    return () => { active = false }
  }, [])
  const selected = useMemo(() => cellTypes.find((cell) => cell.id === selectedId), [cellTypes, selectedId])
  const shownMethods = selected?.methods ?? []
  return (
    <section className="content-page clock-page">
      <p className="eyebrow">Transcriptomic clocks</p>
      <h1>Predicting age from gene expression</h1>
      <p className="clock-intro">Compare a basic expression clock with the most advanced clock available for each cell type. Training performance uses leave-one-study-out predictions; held-out performance comes from studies excluded from model fitting and tuning. The “Across cell types” option averages predictions from separate cell-type clocks for the same donor.</p>
      {status === 'loading' && <p role="status">Loading clock results…</p>}
      {status === 'unavailable' && <p role="status">Clock results are not available in this site release.</p>}
      {status === 'invalid' && <p role="alert">Clock results could not be read because the release data are invalid.</p>}
      {status === 'ready' && cellTypes.length === 0 && <p>No cell types have complete clock results in this release.</p>}
      {cellTypes.length > 0 && (
        <>
          <label className="clock-type-select">Cell type
            <select value={selectedId} onChange={(event) => setSelectedId(event.target.value)}>
              {cellTypes.map((cell) => <option key={cell.id} value={cell.id}>{cell.label}</option>)}
            </select>
          </label>
          <p className="clock-explanation">Each dot is a donor pseudobulk. The diagonal marks perfect agreement between observed and predicted age. Study colors make differences in cohort performance visible.</p>
          {shownMethods.map((method) => <ClockCard key={method.id} method={method} />)}
        </>
      )}
    </section>
  )
}
