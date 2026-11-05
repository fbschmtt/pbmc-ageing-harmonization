import './atlasWorkflow.css'

function WorkflowArt() {
  return (
    <svg viewBox="0 0 720 170" aria-hidden="true" className="workflow-art">
      <defs>
        <marker id="workflow-arrow" viewBox="0 0 8 8" refX="7" refY="4" markerWidth="6" markerHeight="6" orient="auto">
          <path d="M 1 1 L 7 4 L 1 7" fill="none" stroke="#759b8e" strokeWidth="1.4" />
        </marker>
      </defs>
      <path d="M 132 85 H 182 M 338 52 H 365 M 338 120 H 365 M 500 52 H 552 M 500 120 H 552" className="workflow-connector" />
      <path d="M 313 85 H 338 M 338 52 V 120" className="workflow-branch" />
      <g className="workflow-sample" transform="translate(32 27)">
        <rect width="32" height="105" rx="16" />
        <circle cx="12" cy="28" r="5" /><circle cx="21" cy="48" r="5" /><circle cx="12" cy="72" r="5" />
      </g>
      <g className="workflow-sample" transform="translate(72 42)">
        <rect width="32" height="90" rx="16" />
        <circle cx="20" cy="25" r="5" /><circle cx="11" cy="44" r="5" /><circle cx="20" cy="65" r="5" />
      </g>
      <g transform="translate(205 36)">
        <circle cx="25" cy="24" r="19" fill="#e8b196" /><circle cx="79" cy="24" r="19" fill="#a8c6b4" />
        <circle cx="25" cy="78" r="19" fill="#d4c18a" /><circle cx="79" cy="78" r="19" fill="#a6bdce" />
        <g fill="none" stroke="#fffefa" strokeWidth="2" opacity=".8">
          <circle cx="25" cy="24" r="8" /><circle cx="79" cy="24" r="8" />
          <circle cx="25" cy="78" r="8" /><circle cx="79" cy="78" r="8" />
        </g>
      </g>
      <g transform="translate(385 30)" className="workflow-blocks">
        <rect x="0" y="0" width="18" height="18" fill="#e8b196" /><rect x="23" y="0" width="18" height="18" fill="#a8c6b4" /><rect x="46" y="0" width="18" height="18" fill="#d4c18a" /><rect x="69" y="0" width="18" height="18" fill="#a6bdce" />
        <rect x="0" y="24" width="18" height="18" fill="#e8b196" /><rect x="23" y="24" width="18" height="18" fill="#a8c6b4" /><rect x="46" y="24" width="18" height="18" fill="#d4c18a" /><rect x="69" y="24" width="18" height="18" fill="#a6bdce" />
      </g>
      <text x="428" y="20" textAnchor="middle" className="workflow-art-label">GENE COUNTS</text>
      <g transform="translate(385 106)">
        <rect width="87" height="26" rx="13" fill="#e6eee5" />
        <path d="M 0 13 A 13 13 0 0 1 13 0 H 45 V 26 H 13 A 13 13 0 0 1 0 13" fill="#e8b196" />
        <path d="M 45 0 H 67 V 26 H 45 Z" fill="#a8c6b4" />
      </g>
      <text x="428" y="98" textAnchor="middle" className="workflow-art-label">CELL COUNTS</text>
      <g transform="translate(575 25)" fill="none" strokeLinecap="round" strokeLinejoin="round">
        <path d="M 0 52 V 5 M 0 52 H 90" stroke="#bacbc0" />
        <path d="M 5 40 L 22 32 L 39 36 L 56 18 L 76 14 L 88 7" stroke="#bd765e" strokeWidth="3" />
        <path d="M 5 35 L 22 38 L 39 29 L 56 32 L 76 22 L 88 24" stroke="#669789" strokeWidth="3" />
      </g>
      <g transform="translate(575 103)">
        <rect x="0" y="0" width="89" height="11" rx="5.5" fill="#e8b196" />
        <rect x="0" y="17" width="60" height="11" rx="5.5" fill="#a8c6b4" />
        <rect x="0" y="34" width="36" height="11" rx="5.5" fill="#d4c18a" />
      </g>
    </svg>
  )
}

export function AtlasWorkflow() {
  return (
    <section className="atlas-workflow" aria-labelledby="workflow-heading">
      <div className="workflow-intro">
        <div>
          <p className="eyebrow">From studies to insights</p>
          <h2 id="workflow-heading">How the data come together</h2>
        </div>
        <p>Follow the data from individual blood samples to the views below.</p>
      </div>
      <WorkflowArt />
      <ol className="workflow-steps">
        <li><span>01</span><strong>Collect studies</strong><p>PBMC data from individual samples and studies.</p></li>
        <li><span>02</span><strong>Harmonize cell types</strong><p>Assign a shared AIFI annotation to cells.</p></li>
        <li><span>03</span><strong>Summarize samples</strong><p>Sum gene counts by sample and cell type; count cells for composition.</p></li>
        <li>
          <span>04</span><strong>Explore patterns</strong>
          <p>Compare results across age, studies, and cell types.</p>
          <div className="workflow-links" aria-label="Explore atlas results">
            <a href="#differential">Expression</a>
            <a href="#composition">Composition</a>
            <a href="#trajectories">Trajectories</a>
          </div>
        </li>
      </ol>
      <p className="workflow-note">Expression models use sample-level pseudobulks. Composition describes relative cell fractions.</p>
    </section>
  )
}
