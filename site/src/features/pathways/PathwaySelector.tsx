import { GeneSetInput } from '../../components/GeneSetInput'
import { PathwayAttribution } from '../../components/PathwayAttribution'
import { parseGeneSet } from '../../data/geneSet'
import type { PathwaySet } from '../../data/pathways'
import './pathways.css'

type Props = {
  genes: string
  onGenesChange: (genes: string) => void
  pathwaySets: PathwaySet[]
  activePathway: PathwaySet | null
  onSelectPathway: (pathway: PathwaySet) => void
  onOpenMatrix: () => void
}

export function PathwaySelector({
  genes,
  onGenesChange,
  pathwaySets,
  activePathway,
  onSelectPathway,
  onOpenMatrix,
}: Props) {
  const count = parseGeneSet(genes).size
  return (
    <section className="pathways-page pathway-selector">
      <h1 className="visually-hidden">Choose genes to highlight</h1>
      <div className="pathway-selector-card">
        <GeneSetInput
          label="Gene set"
          value={genes}
          onChange={onGenesChange}
          helper="Enter symbols separated by commas or whitespace. Matching is exact and case-insensitive."
          pathwaySets={pathwaySets}
          activePathway={activePathway}
          onSelectPathway={onSelectPathway}
          pathwayDisplay="inline"
          onClear={() => onGenesChange('')}
        />
        {pathwaySets.length > 0 && <PathwayAttribution activePathway={activePathway} />}
        <p className="pathway-count">
          {count.toLocaleString()} {count === 1 ? 'gene' : 'genes'} selected
        </p>
      </div>
      <p className="pathway-intro">
        The selection highlights matching genes in differential expression, cell-type profiles and the
        cross-type trajectory map. It also defines rows in the age-bin gene × cell type matrix. It does not
        change the fitted analyses or perform pathway enrichment.
      </p>
      <button className="primary" onClick={onOpenMatrix}>
        Open age-bin gene × cell type →
      </button>
    </section>
  )
}
