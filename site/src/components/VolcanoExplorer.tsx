import type { ReactNode } from 'react'
import type { PathwaySet } from '../data/pathways'
import { GeneSetInput } from './GeneSetInput'
import './volcanoExplorer.css'

type GeneSetControl = {
  value: string
  onChange: (value: string) => void
  pathwaySets: PathwaySet[]
  activePathway: PathwaySet | null
  onSelectPathway?: (pathway: PathwaySet) => void
}

export function VolcanoExplorer({
  controls,
  geneSet,
  title,
  description,
  legend,
  plot,
  detail,
  descriptionPlacement = 'top',
  legendPlacement = 'top',
}: {
  controls: ReactNode
  geneSet?: GeneSetControl
  title: string
  description: ReactNode
  legend: ReactNode
  plot: ReactNode
  detail: ReactNode
  descriptionPlacement?: 'top' | 'below'
  legendPlacement?: 'top' | 'below'
}) {
  return (
    <div className="volcano-explorer">
      <aside className="volcano-explorer-controls">
        {controls}
        {geneSet && (
          <details className="volcano-explorer-gene-set" open>
            <summary>
              Gene sets{geneSet.value.trim() ? <span> · active</span> : null}
            </summary>
            <GeneSetInput
              label="Highlight genes"
              value={geneSet.value}
              onChange={geneSet.onChange}
              onClear={() => geneSet.onChange('')}
              helper="This gene set is shared across the atlas."
              pathwaySets={geneSet.pathwaySets}
              activePathway={geneSet.activePathway}
              onSelectPathway={geneSet.onSelectPathway}
            />
          </details>
        )}
      </aside>
      <div className="volcano-explorer-main">
        <p className="eyebrow">Visualise</p>
        <h2>{title}</h2>
        {descriptionPlacement === 'top' && (
          <div className="volcano-explorer-description">{description}</div>
        )}
        {legendPlacement === 'top' && <div className="volcano-explorer-legend">{legend}</div>}
        {plot}
        {legendPlacement === 'below' && <div className="volcano-explorer-legend">{legend}</div>}
        {descriptionPlacement === 'below' && (
          <details className="volcano-explorer-description volcano-explorer-about">
            <summary>About this plot</summary>
            {description}
          </details>
        )}
      </div>
      <aside className="volcano-explorer-detail">{detail}</aside>
    </div>
  )
}
