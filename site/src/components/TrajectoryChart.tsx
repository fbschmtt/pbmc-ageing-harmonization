import type { TrajectoryRecord } from '../data/atlas'

export function TrajectoryChart({ record }: { record: TrajectoryRecord }) {
  const values = record.trajectory
  const extent = Math.max(0.25, ...values.map(Math.abs))
  const plotStart = 14
  const plotEnd = 272
  const pointSpacing = (plotEnd - plotStart) / (values.length - 1)
  const points = values
    .map((value, index) => `${plotStart + index * pointSpacing},${48 - (value / extent) * 36}`)
    .join(' ')

  return (
    <svg
      className="trajectory-chart"
      viewBox="0 0 286 96"
      role="img"
      aria-label={`Scaled age-bin effect trajectory for ${record.gene} in ${record.cellType}`}
    >
      <line x1={plotStart} x2={plotEnd} y1="48" y2="48" className="chart-zero" strokeDasharray="3 3" />
      <polyline
        points={points}
        className={values.at(-1)! >= 0 ? 'trajectory-line rise' : 'trajectory-line fall'}
      />
      {values.map((value, index) => (
        <circle
          key={index}
          cx={plotStart + index * pointSpacing}
          cy={48 - (value / extent) * 36}
          r="2.5"
          className="trajectory-dot"
        />
      ))}
    </svg>
  )
}
