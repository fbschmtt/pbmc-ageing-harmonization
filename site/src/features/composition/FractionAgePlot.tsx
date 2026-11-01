export type { FractionAgeEstimate } from '../../data/contracts.mjs'
import type { FractionAgeEstimate } from '../../data/contracts.mjs'

export function FractionAgePlot({
  estimates,
  studyColors,
}: {
  estimates: FractionAgeEstimate[]
  studyColors: ReadonlyMap<string, string>
}) {
  if (!estimates.length) return <p>No estimable age-bin fraction model is available for this cell type.</p>
  const groups = [...new Set(estimates.map((item) => item.group))].sort()
  const bins = [...new Set(estimates.map((item) => item.ageBin))].sort(
    (a, b) => Number(a.split('-')[0]) - Number(b.split('-')[0]),
  )
  const points = estimates.filter((item) => [item.estimate, item.low, item.high].every(Number.isFinite))
  if (!points.length) return <p>No finite age-bin estimates are available for this cell type.</p>
  const minimum = Math.min(0, ...points.map((item) => item.low))
  const maximum = Math.max(0, ...points.map((item) => item.high))
  const padding = Math.max((maximum - minimum) * 0.1, 0.01)
  const low = minimum - padding,
    high = maximum + padding
  const left = 70,
    right = 790,
    top = 25,
    bottom = 300
  const x = (index: number) => left + (index * (right - left)) / Math.max(1, bins.length - 1)
  const y = (value: number) => bottom - ((value - low) / (high - low)) * (bottom - top)
  return (
    <div className="fraction-chart-panel">
      <p className="eyebrow">Age profile</p>
      <h2>Adjusted fraction change across age</h2>
      <p>
        Each study is aligned to zero in its {estimates[0].referenceBin} reference group. Points show adjusted
        percentage-point differences, with pointwise 95% confidence intervals and sample counts.
      </p>
      <div className="fraction-chart-scroll">
        <svg
          viewBox="0 0 820 360"
          role="img"
          aria-label="Adjusted fraction changes by age bin with 95 percent confidence intervals"
        >
          {[0, 0.25, 0.5, 0.75, 1].map((step) => {
            const value = low + step * (high - low)
            return (
              <g key={step}>
                <line x1={left} x2={right} y1={y(value)} y2={y(value)} stroke="#e1e7df" />
                <text x={left - 10} y={y(value) + 4} textAnchor="end" className="fraction-chart-text">
                  {(value * 100).toFixed(1)} pp
                </text>
              </g>
            )
          })}
          {bins.map((bin, index) => (
            <text key={bin} x={x(index)} y={bottom + 25} textAnchor="middle" className="fraction-chart-text">
              {bin}
            </text>
          ))}
          {groups.flatMap((group, groupIndex) => {
            const offset = (groupIndex - (groups.length - 1) / 2) * Math.min(10, 30 / groups.length)
            const byBin = new Map(
              points.filter((item) => item.group === group).map((item) => [item.ageBin, item]),
            )
            return bins.slice(1).flatMap((bin, index) => {
              const before = byBin.get(bins[index])
              const after = byBin.get(bin)
              const beforeStart = Number(bins[index].split('-')[0])
              const beforeEnd = Number(bins[index].split('-')[1])
              if (!before || !after || Number(bin.split('-')[0]) !== beforeEnd || beforeEnd <= beforeStart)
                return []
              return [
                <line
                  key={`${group}-${bin}`}
                  x1={x(index) + offset}
                  x2={x(index + 1) + offset}
                  y1={y(before.estimate)}
                  y2={y(after.estimate)}
                  stroke={studyColors.get(group) ?? '#65736f'}
                  strokeWidth="1.5"
                />,
              ]
            })
          })}
          {points.map((item, index) => {
            const group = groups.indexOf(item.group),
              offset = (group - (groups.length - 1) / 2) * Math.min(10, 30 / groups.length),
              px = x(bins.indexOf(item.ageBin)) + offset,
              color = studyColors.get(item.group) ?? '#65736f'
            return (
              <g key={`${item.group}-${item.ageBin}-${index}`}>
                <line x1={px} x2={px} y1={y(item.low)} y2={y(item.high)} stroke={color} strokeWidth="2" />
                <circle cx={px} cy={y(item.estimate)} r="4" fill={color}>
                  <title>
                    {item.group}, {item.ageBin}: {(item.estimate * 100).toFixed(2)} pp (95% CI{' '}
                    {(item.low * 100).toFixed(2)}–{(item.high * 100).toFixed(2)} pp), n={item.samples}
                  </title>
                </circle>
              </g>
            )
          })}
          <text x={(left + right) / 2} y="350" textAnchor="middle" className="fraction-chart-text">
            Age (years)
          </text>
        </svg>
      </div>
      <div className="fraction-legend">
        {groups.map((group) => (
          <span key={group}>
            <i style={{ backgroundColor: studyColors.get(group) ?? '#65736f' }} />
            {group}
          </span>
        ))}
      </div>
    </div>
  )
}
