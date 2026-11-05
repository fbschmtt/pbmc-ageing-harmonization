import type { Page } from '../../app/navigation'
import { AtlasWorkflow } from './AtlasWorkflow'

export function OverviewPage({
  onNavigate,
  summary,
}: {
  onNavigate: (page: Page) => void
  summary?: { studies: number; donors: number; cells: number } | null
}) {
  return (
    <>
      <section className="hero">
        <p className="eyebrow">Human peripheral blood</p>
        <h1>Ageing across immune cell types</h1>
        <p className="hero-copy">
          This atlas brings together single-cell studies of peripheral blood mononuclear cells (PBMCs) to
          examine how gene expression and cell-type composition vary with age.
        </p>
        {summary && (
          <dl className="atlas-summary" aria-label="Single-cell atlas coverage">
            <div>
              <dt>Studies</dt>
              <dd>{summary.studies.toLocaleString()}</dd>
            </div>
            <div>
              <dt>Individual donors</dt>
              <dd>{summary.donors.toLocaleString()}</dd>
            </div>
            <div>
              <dt>Cells</dt>
              <dd>{summary.cells.toLocaleString()}</dd>
            </div>
          </dl>
        )}
        <p className="hero-copy">
          Cells are assigned a shared cell-type annotation, then aggregated by sample and cell type for
          analysis. Differences between studies and uneven sample coverage should be considered when
          interpreting the results.
        </p>
        <div className="actions">
          <button className="primary" onClick={() => onNavigate('differential')}>
            View differential expression <span>→</span>
          </button>
          <button className="secondary" onClick={() => onNavigate('modelling')}>
            About the analysis
          </button>
        </div>
      </section>
      <AtlasWorkflow />
      <section className="entry-grid">
        <article>
          <a href="#differential">
            <h2>Differential expression</h2>
            <p>Examine associations between gene expression and age or other covariates within cell types.</p>
            <span>View expression results →</span>
          </a>
        </article>
        <article>
          <a href="#composition">
            <h2>Cell composition</h2>
            <p>Compare how cell-type fractions vary across studies and age groups.</p>
            <span>View composition →</span>
          </a>
        </article>
        <article>
          <a href="#trajectories">
            <h2>Age-bin trajectories</h2>
            <p>Compare patterns across age groups for genes selected for trajectory clustering.</p>
            <span>View trajectory map →</span>
          </a>
        </article>
      </section>
    </>
  )
}
