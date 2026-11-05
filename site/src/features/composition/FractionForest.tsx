export type { FractionCoefficient } from '../../data/contracts.mjs'
import type { FractionCoefficient } from '../../data/contracts.mjs'
import { formatAxisValue, symmetricAxis } from '../../theme/axisScale'

export type FractionDiagnostic = {
  study: string
  residualSd: number | null
  binomialSamplingSd: number | null
}

const names: Record<string, string> = {
  age: 'Age',
  sex: 'Sex',
  bmi: 'BMI',
  cmv: 'CMV',
  smoking_status: 'Smoking status',
}
const order = ['age', 'sex', 'bmi', 'cmv', 'smoking_status']

export function FractionForest({
  coefficients,
  diagnostics,
  studyColors,
}: {
  coefficients: FractionCoefficient[]
  diagnostics: FractionDiagnostic[]
  studyColors: ReadonlyMap<string, string>
}) {
  if (!coefficients.length)
    return (
        <p>No coefficient estimates are available for this cell type.</p>
    )
  const studies = [...new Set([...coefficients, ...diagnostics].map((item) => item.study))].sort()
  const hasDiagnostics = diagnostics.some(
    (item) => item.residualSd != null || item.binomialSamplingSd != null,
  )
  const rows = [...new Set(coefficients.map((item) => `${item.covariate}\0${item.contrast}`))].sort(
    (a, b) => {
      const [aName] = a.split('\0'),
        [bName] = b.split('\0')
      return (
        (order.indexOf(aName) < 0 ? 99 : order.indexOf(aName)) -
          (order.indexOf(bName) < 0 ? 99 : order.indexOf(bName)) || a.localeCompare(b)
      )
    },
  )
  const axis = symmetricAxis([
      0.01,
      ...coefficients.flatMap((item) => [item.low, item.high, item.estimate]),
      ...diagnostics.flatMap((item) => [item.residualSd ?? 0, item.binomialSamplingSd ?? 0]),
    ], 0.1)
  const extent = axis.extent
  const width = 820,
    left = 190,
    right = 790,
    top = 45,
    height = Math.max(190, (rows.length + (hasDiagnostics ? 1 : 0)) * 65 + 110)
  const x = (value: number) => left + ((value / extent + 1) / 2) * (right - left)
  const rowY = (index: number) => top + index * 65 + 20
  return (
    <div className="fraction-chart-panel atlas-card atlas-card-soft">
      <div className="fraction-chart-heading">
        <div>
          <p className="eyebrow">Study-specific models</p>
          <h2>Covariate effects</h2>
          <p>
            Adjusted effects on fraction within {coefficients[0].parent}. Points show estimates; lines show
            95% confidence intervals. Residual SD markers compare the fitted model's residual variation with
            the conditional binomial sampling reference; they are not confidence intervals.
          </p>
        </div>
      </div>
      <div className="fraction-chart-scroll">
        <svg
          viewBox={`0 0 ${width} ${height}`}
          role="img"
          aria-label="Forest plot of study-specific covariate effects on parent-compartment fractions"
        >
          <line x1={x(0)} x2={x(0)} y1={top - 10} y2={height - 55} stroke="#8fa39a" strokeDasharray="4 5" />
          {rows.map((key, index) => {
            const [covariate, contrast] = key.split('\0')
            return (
              <g key={key}>
                <text x={left - 15} y={rowY(index) + 4} textAnchor="end" className="fraction-chart-text">
                  {names[covariate] || covariate}: {contrast}
                </text>
                <line x1={left} x2={right} y1={rowY(index)} y2={rowY(index)} stroke="#edf0ea" />
              </g>
            )
          })}
          {hasDiagnostics && (
            <g>
              <text x={left - 15} y={rowY(rows.length) + 4} textAnchor="end" className="fraction-chart-text">
                Residual SD
              </text>
              <line x1={left} x2={right} y1={rowY(rows.length)} y2={rowY(rows.length)} stroke="#edf0ea" />
            </g>
          )}
          {coefficients.map((item, index) => {
            const row = rows.indexOf(`${item.covariate}\0${item.contrast}`),
              study = studies.indexOf(item.study),
              offset = (study - (studies.length - 1) / 2) * Math.min(8, 34 / Math.max(studies.length, 1)),
              cy = rowY(row) + offset,
              color = studyColors.get(item.study) ?? '#65736f'
            return (
              <g key={`${item.study}-${item.term}-${index}`}>
                <line x1={x(item.low)} x2={x(item.high)} y1={cy} y2={cy} stroke={color} strokeWidth="2" />
                <line x1={x(item.low)} x2={x(item.low)} y1={cy - 3} y2={cy + 3} stroke={color} />
                <line x1={x(item.high)} x2={x(item.high)} y1={cy - 3} y2={cy + 3} stroke={color} />
                <circle cx={x(item.estimate)} cy={cy} r="4" fill={color}>
                  <title>
                    {item.study}: {(item.estimate * 100).toFixed(2)} pp (95% CI {(item.low * 100).toFixed(2)}{' '}
                    to {(item.high * 100).toFixed(2)}); n={item.samples}
                  </title>
                </circle>
              </g>
            )
          })}
          {diagnostics.map((item) => {
            const study = studies.indexOf(item.study)
            const offset = (study - (studies.length - 1) / 2) * Math.min(8, 34 / Math.max(studies.length, 1))
            const cy = rowY(rows.length) + offset
            const color = studyColors.get(item.study) ?? '#65736f'
            return (
              <g key={`diagnostic-${item.study}`}>
                {item.residualSd != null && (
                  <circle cx={x(item.residualSd)} cy={cy} r="4" fill={color}>
                    <title>
                      {item.study}: observed residual SD {(item.residualSd * 100).toFixed(3)} pp
                    </title>
                  </circle>
                )}
                {item.binomialSamplingSd != null && (
                  <path
                    d={`M ${x(item.binomialSamplingSd)} ${cy - 5} l 5 5 l -5 5 l -5 -5 z`}
                    fill="none"
                    stroke={color}
                    strokeWidth="1.5"
                  >
                    <title>
                      {item.study}: simulated sampling SD {(item.binomialSamplingSd * 100).toFixed(3)} pp
                    </title>
                  </path>
                )}
              </g>
            )
          })}
          {axis.ticks.map((tick) => (
            <text
              key={tick}
              x={x(tick)}
              y={height - 32}
              textAnchor="middle"
              className="fraction-chart-text"
              fontWeight={tick === 0 ? '700' : undefined}
            >
              {formatAxisValue(tick * 100, extent * 100)}
            </text>
          ))}
          <text x={(left + right) / 2} y={height - 8} textAnchor="middle" className="fraction-chart-text">
            Fraction change (percentage points)
          </text>
        </svg>
      </div>
      <div className="fraction-legend">
        {studies.map((study) => (
          <span key={study}>
            <i style={{ backgroundColor: studyColors.get(study) ?? '#65736f' }} />
            {study}
          </span>
        ))}
      </div>
      {hasDiagnostics && (
        <p className="fraction-diagnostic-key">● Observed residual SD &nbsp;&nbsp; ◇ Simulated sampling SD</p>
      )}
    </div>
  )
}
