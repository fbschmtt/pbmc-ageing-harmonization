import { loadAgePreview, type Bin, type Preview } from '../../data/agePreview'
import { visualColors } from '../../theme/visualEncoding'
import { useEffect, useState } from 'react'

export function GeneAgePreview({
  gene,
  cellType,
  asset,
}: {
  gene: string
  cellType: string
  asset: string | null
}) {
  const [loaded, setLoaded] = useState<{ asset: string; gene: string; data: Preview } | null>(null)
  const [failed, setFailed] = useState(false)
  useEffect(() => {
    if (!asset) return
    let active = true
    setFailed(false)
    loadAgePreview(asset, gene)
      .then((data) => {
        if (active) setLoaded({ asset, gene, data })
      })
      .catch(() => {
        if (active) setFailed(true)
      })
    return () => {
      active = false
    }
  }, [asset, gene])
  if (!asset) return <p>Combined age-bin model is unavailable for this cell type.</p>
  if (failed) return <p>Age-bin results could not be loaded.</p>
  if (loaded?.asset !== asset || loaded.gene !== gene) return <p>Loading age-bin results…</p>
  const bins = loaded.data.genes[gene]
  if (!Array.isArray(bins) || !bins.length)
    return <p>No combined age-bin result is available for this gene.</p>
  const valid = bins.filter(
    (item): item is Bin & { effect: number } => item.effect !== null && Number.isFinite(item.effect),
  )
  if (!valid.length) return <p>Age-bin estimates are unavailable for this gene.</p>
  const extent = Math.max(0.1, ...valid.map((item) => Math.abs(item.effect))) * 1.15
  const width = 300,
    left = 38,
    right = 290,
    top = 15,
    bottom = 125
  const starts = bins.map((item) => Number(item.bin.split('-')[0]))
  const first = Math.min(...starts),
    last = Math.max(...starts)
  const x = (start: number) => left + ((start - first) / Math.max(1, last - first)) * (right - left)
  const y = (value: number) => top + ((extent - value) / (2 * extent)) * (bottom - top)
  return (
    <div className="de-trajectory-preview">
      <h3>Expression across age</h3>
      <p>
        Age-bin fit relative to {loaded.data.referenceBin || 'the reference bin'}. This is separate from the
        continuous-age analysis selected in the plot.
      </p>
      <svg
        viewBox={`0 0 ${width} 160`}
        role="img"
        aria-label={`Combined age-bin effects for ${gene} in ${cellType}, relative to ${loaded.data.referenceBin || 'the reference bin'}`}
      >
        <line x1={left} x2={right} y1={y(0)} y2={y(0)} stroke="#aebeb3" strokeDasharray="3 3" />
        {[-extent, 0, extent].map((tick) => (
          <text key={tick} x={left - 6} y={y(tick) + 3} textAnchor="end">
            {tick.toFixed(1)}
          </text>
        ))}
        {bins.slice(1).map((current, index) => {
          const previous = bins[index]
          if (
            previous.effect === null ||
            current.effect === null ||
            !Number.isFinite(previous.effect) ||
            !Number.isFinite(current.effect) ||
            Number(previous.bin.split('-')[1]) !== Number(current.bin.split('-')[0])
          )
            return null
          return (
            <line
              key={current.bin}
              x1={x(starts[index])}
              x2={x(starts[index + 1])}
              y1={y(previous.effect)}
              y2={y(current.effect)}
              stroke={visualColors.brain}
              strokeWidth="2"
            />
          )
        })}
        {valid.map((item) => (
          <circle
            key={item.bin}
            cx={x(Number(item.bin.split('-')[0]))}
            cy={y(item.effect)}
            r="3"
            fill={visualColors.brain}
          >
            <title>
              {item.bin}: {item.effect.toFixed(3)} log₂ fold change
            </title>
          </circle>
        ))}
        {bins.map((item) => (
          <text key={item.bin} x={x(Number(item.bin.split('-')[0]))} y="143" textAnchor="middle">
            {item.bin}
          </text>
        ))}
      </svg>
      {loaded.data.design && <p>Age-bin model: {loaded.data.design}</p>}
    </div>
  )
}
