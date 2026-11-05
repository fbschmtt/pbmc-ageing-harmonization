export type { FractionAgeEstimate } from '../../data/contracts.mjs'
import type { FractionAgeEstimate } from '../../data/contracts.mjs'
import { formatAxisValue } from '../../theme/axisScale'
import { absoluteFractionPoints, clipFractionInterval, pointBasedFractionAxis } from './ageFractionPlot'

export function FractionAgePlot({
  estimates,
  studyColors,
}: {
  estimates: FractionAgeEstimate[]
  studyColors: ReadonlyMap<string, string>
}) {
  if (!estimates.length) return <p>No estimable age-bin fraction model is available for this cell type.</p>
  const points = absoluteFractionPoints(estimates)
  if (!points.length) return <p>Absolute reference-bin estimates are unavailable.</p>
  const groups = [...new Set(points.map(({ item }) => item.group))].sort()
  const bins = [...new Set(estimates.map((item) => item.ageBin))].sort(
    (a, b) => Number(a.split('-')[0]) - Number(b.split('-')[0]),
  )
  const { low, high, ticks } = pointBasedFractionAxis(points)
  const left = 70,
    right = 790,
    top = 25,
    bottom = 300
  const x = (index: number) => left + (index * (right - left)) / Math.max(1, bins.length - 1)
  const y = (value: number) => bottom - ((value - low) / (high - low)) * (bottom - top)
  return (
    <div className="fraction-chart-panel atlas-card atlas-card-soft">
      <p className="eyebrow">Age profile</p>
      <h2>Adjusted fraction across age</h2>
      <p>
        Points show estimated absolute fractions, anchored at each study's model-adjusted {estimates[0].referenceBin}
        {' '}reference fraction. Other bins add their estimated age-bin contrast to that reference level. Whiskers show
        pointwise 95% contrast intervals shifted to the reference level; intervals extending beyond the
        point-based y-axis are clipped and marked with arrows.
      </p>
      <div className="fraction-chart-scroll">
        <svg
          viewBox="0 0 820 360"
          role="img"
          aria-label="Model-adjusted absolute fractions by age bin with clipped 95 percent confidence intervals"
        >
          {ticks.map((value) => {
            return (
              <g key={value}>
                <line x1={left} x2={right} y1={y(value)} y2={y(value)} stroke={value === 0 ? '#8fa39a' : '#e1e7df'} strokeDasharray={value === 0 ? '4 5' : undefined} />
                <text x={left - 10} y={y(value) + 4} textAnchor="end" className="fraction-chart-text" fontWeight={value === 0 ? '700' : undefined}>
                  {formatAxisValue(value * 100, Math.max(Math.abs(low), Math.abs(high)) * 100)}%
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
              points.filter(({ item }) => item.group === group).map(({ item, ...point }) => [item.ageBin, { item, ...point }]),
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
          {points.map(({ item, estimate, low: intervalLow, high: intervalHigh }, index) => {
            const group = groups.indexOf(item.group),
              offset = (group - (groups.length - 1) / 2) * Math.min(10, 30 / groups.length),
              px = x(bins.indexOf(item.ageBin)) + offset,
              color = studyColors.get(item.group) ?? '#65736f',
              interval = clipFractionInterval({ item, estimate, low: intervalLow, high: intervalHigh }, low, high)
            return (
              <g key={`${item.group}-${item.ageBin}-${index}`}>
                <line x1={px} x2={px} y1={y(interval.low)} y2={y(interval.high)} stroke={color} strokeWidth="2" />
                {interval.clippedLow && (
                  <path d={`M ${px - 4} ${y(low) - 5} L ${px + 4} ${y(low) - 5} L ${px} ${y(low) + 1} Z`} fill={color} />
                )}
                {interval.clippedHigh && (
                  <path d={`M ${px - 4} ${y(high) + 5} L ${px + 4} ${y(high) + 5} L ${px} ${y(high) - 1} Z`} fill={color} />
                )}
                <circle cx={px} cy={y(estimate)} r="4" fill={color}>
                  <title>
                    {item.group}, {item.ageBin}: {(estimate * 100).toFixed(2)}% estimated fraction (shifted 95% CI{' '}
                    {(intervalLow * 100).toFixed(2)}–{(intervalHigh * 100).toFixed(2)}%), n={item.samples}
                  </title>
                </circle>
              </g>
            )
          })}
          <text x={(left + right) / 2} y="350" textAnchor="middle" className="fraction-chart-text">
            Age (years)
          </text>
          <text
            x="15"
            y={(top + bottom) / 2}
            textAnchor="middle"
            className="fraction-chart-text"
            transform={`rotate(-90 15 ${(top + bottom) / 2})`}
          >
            Adjusted fraction (%)
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
