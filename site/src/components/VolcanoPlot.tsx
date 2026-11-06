import { formatAxisValue, positiveAxis, symmetricAxis } from '../theme/axisScale'
import './volcanoPlot.css'

type Appearance = { fill: string; stroke?: string; strokeWidth?: number; radius?: number; opacity?: number }
type Props<T> = {
  records: T[]
  selected: T | null
  effect: (record: T) => number
  adjustedP: (record: T) => number
  id: (record: T) => string
  title: (record: T) => string
  appearance: (record: T, selected: boolean) => Appearance
  highlighted?: (record: T) => boolean
  onSelect: (record: T) => void
  xLabel: string
  ariaLabel: string
  threshold?: number
  className?: string
}

const negLog10 = (value: number) => -Math.log10(Math.max(value, Number.MIN_VALUE))
const width = 800,
  height = 600,
  left = 72,
  right = 772,
  top = 26,
  bottom = 542

export function VolcanoPlot<T>({
  records,
  selected,
  effect,
  adjustedP,
  id,
  title,
  appearance,
  highlighted,
  onSelect,
  xLabel,
  ariaLabel,
  threshold,
  className = '',
}: Props<T>) {
  if (!records.length)
    return <div className="empty-scatter">No plottable results are available for this selection.</div>
  const xAxis = symmetricAxis(records.map(effect))
  const yAxis = positiveAxis([
    ...records.map((record) => negLog10(adjustedP(record))),
    ...(threshold ? [negLog10(threshold)] : []),
  ])
  const x = (value: number) => left + ((value + xAxis.extent) / (2 * xAxis.extent)) * (right - left)
  const y = (value: number) => bottom - (value / yAxis.maximum) * (bottom - top)
  const selectedIndex = selected === null ? -1 : records.indexOf(selected)
  const ordered = records.map((record, index) => ({ record, index })).sort(
    (a, b) =>
      Number(a.index === selectedIndex) - Number(b.index === selectedIndex) ||
      Number(Boolean(highlighted?.(a.record))) - Number(Boolean(highlighted?.(b.record))),
  )
  return (
    <figure className={`volcano-plot ${className}`}>
      <svg viewBox={`0 0 ${width} ${height}`} role="img" aria-label={ariaLabel}>
        {yAxis.ticks.map((tick) => (
          <g key={tick}>
            <line x1={left} y1={y(tick)} x2={right} y2={y(tick)} className="volcano-grid" />
            <text x={left - 12} y={y(tick) + 4} textAnchor="end">
              {formatAxisValue(tick, yAxis.maximum)}
            </text>
          </g>
        ))}
        <line x1={x(0)} y1={top} x2={x(0)} y2={bottom} className="volcano-zero" />
        {threshold && (
          <line
            x1={left}
            y1={y(negLog10(threshold))}
            x2={right}
            y2={y(negLog10(threshold))}
            className="volcano-threshold"
          />
        )}
        {xAxis.ticks.map((tick) => (
          <text key={tick} x={x(tick)} y={bottom + 22} textAnchor="middle">
            {formatAxisValue(tick, xAxis.extent)}
          </text>
        ))}
        <text x={(left + right) / 2} y={height - 10} textAnchor="middle">
          {xLabel}
        </text>
        <text
          x="18"
          y={(top + bottom) / 2}
          textAnchor="middle"
          transform={`rotate(-90 18 ${(top + bottom) / 2})`}
        >
          −log₁₀ adjusted p-value
        </text>
        {ordered.map(({ record, index }) => {
          const isSelected = index === selectedIndex
          const style = appearance(record, isSelected)
          return (
            <circle
              key={`${id(record)}-${index}`}
              cx={x(effect(record))}
              cy={y(negLog10(adjustedP(record)))}
              r={style.radius ?? 5}
              fill={style.fill}
              stroke={style.stroke ?? 'none'}
              strokeWidth={style.strokeWidth ?? 0}
              opacity={style.opacity ?? 1}
              onClick={() => onSelect(record)}
            >
              <title>{title(record)}</title>
            </circle>
          )
        })}
      </svg>
    </figure>
  )
}
