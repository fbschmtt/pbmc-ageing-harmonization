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
}: {
  controls: ReactNode
  geneSet?: GeneSetControl
  title: string
  description: ReactNode
  legend: ReactNode
  plot: ReactNode
  detail: ReactNode
}) {
  return (
    <div className="volcano-explorer">
      <aside className="volcano-explorer-controls">
        {controls}
        {geneSet && (
          <GeneSetInput
            className="volcano-explorer-gene-set"
            label="Highlight genes"
            value={geneSet.value}
            onChange={geneSet.onChange}
            onClear={() => geneSet.onChange('')}
            helper="This gene set is shared across the atlas."
            pathwaySets={geneSet.pathwaySets}
            activePathway={geneSet.activePathway}
            onSelectPathway={geneSet.onSelectPathway}
          />
        )}
      </aside>
      <div className="volcano-explorer-main">
        <h2>{title}</h2>
        <div className="volcano-explorer-description">{description}</div>
        <div className="volcano-explorer-legend">{legend}</div>
        {plot}
      </div>
      <aside className="volcano-explorer-detail">{detail}</aside>
    </div>
  )
}
