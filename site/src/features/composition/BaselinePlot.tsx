import type { StudyBaseline } from './baselineSummary'

export function BaselinePlot({
  studies,
  mean,
  studyColors,
}: {
  studies: StudyBaseline[]
  mean: number
  studyColors: ReadonlyMap<string, string>
}) {
  const left = 180,
    right = 700,
    top = 22,
    bottom = top + studies.length * 35 + 12
  const minimum = Math.min(0, ...studies.map((item) => item.value), mean)
  const maximum = Math.max(1, ...studies.map((item) => item.value), mean)
  const padding = (maximum - minimum) * 0.05
  const lo = minimum - padding,
    hi = maximum + padding
  const x = (value: number) => left + ((value - lo) / (hi - lo)) * (right - left)
  return (
    <div className="baseline-plot-scroll">
      <svg
        viewBox={`0 0 760 ${bottom + 47}`}
        role="img"
        aria-label="Model-adjusted reference fractions by study, with an equal-weight study mean"
      >
        <line x1={x(mean)} x2={x(mean)} y1={top - 10} y2={bottom} stroke="#496b60" strokeDasharray="4 4">
          <title>Equal-weight mean: {(mean * 100).toFixed(2)}%</title>
        </line>
        {studies.map((item, index) => {
          const y = top + index * 35 + 8
          return (
            <g key={item.study}>
              <text x={left - 12} y={y + 4} textAnchor="end" className="fraction-chart-text">
                {item.study}
              </text>
              <line x1={left} x2={right} y1={y} y2={y} stroke="#e1e7df" />
              <circle cx={x(item.value)} cy={y} r="5" fill={studyColors.get(item.study) ?? '#65736f'}>
                <title>
                  {item.study}: {(item.value * 100).toFixed(2)}%; reference n={item.samples}
                </title>
              </circle>
              <text x={right + 8} y={y + 4} className="fraction-chart-text">
                {(item.value * 100).toFixed(1)}%
              </text>
            </g>
          )
        })}
        {[0, 0.25, 0.5, 0.75, 1].map((tick) => (
          <text key={tick} x={x(tick)} y={bottom + 18} textAnchor="middle" className="fraction-chart-text">
            {(tick * 100).toFixed(0)}%
          </text>
        ))}
        <text x={(left + right) / 2} y={bottom + 38} textAnchor="middle" className="fraction-chart-text">
          Adjusted reference fraction
        </text>
      </svg>
    </div>
  )
}
