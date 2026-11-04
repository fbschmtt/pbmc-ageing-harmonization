import './ageModelGuide.css'

function MiniPlot({ kind }: { kind: 'adult' | 'juvenile' | 'bins' }) {
  return (
    <svg viewBox="0 0 300 150" className="age-model-plot" aria-hidden="true">
      <path d="M 28 15 V 117 H 282" className="age-model-axis" />
      {kind === 'adult' && <>
        <path d="M 80 91 L 258 37" className="age-model-line adult" />
        <path d="M 80 117 V 123 M 258 117 V 123" className="age-model-axis" />
        <text x="80" y="139" textAnchor="middle">20</text><text x="258" y="139" textAnchor="middle">Older</text>
      </>}
      {kind === 'juvenile' && <>
        <path d="M 138 16 V 117" className="age-model-boundary" />
        <path d="M 38 95 L 135 52" className="age-model-line juvenile" />
        <path d="M 140 67 L 267 45" className="age-model-line context" />
        <text x="38" y="139" textAnchor="middle">Younger</text><text x="138" y="139" textAnchor="middle">20</text><text x="266" y="139" textAnchor="middle">Older</text>
      </>}
      {kind === 'bins' && <>
        <path d="M 39 81 H 72 M 73 75 H 106 M 107 82 H 140 M 141 67 H 174 M 175 57 H 208 M 209 64 H 242 M 243 48 H 276" className="age-model-bins" />
        <path d="M 39 81 H 276" className="age-model-reference" />
        <circle cx="55" cy="81" r="4" className="age-model-reference-dot" />
        <text x="55" y="139" textAnchor="middle">20–30</text><text x="259" y="139" textAnchor="middle">80–90</text>
      </>}
      <text x="282" y="12" textAnchor="end">Expression</text>
    </svg>
  )
}

export function AgeModelGuide() {
  return (
    <section className="age-model-guide" aria-labelledby="age-model-heading">
      <div className="age-model-intro">
        <p className="eyebrow">A guide to the age models</p>
        <h2 id="age-model-heading">Three ways to look at age</h2>
        <p>Different models ask different questions of the same kind of sample-level data. These illustrations are examples, not results.</p>
      </div>
      <div className="age-model-grid">
        <article className="age-model-card">
          <span className="age-model-number">01 / THE ADULT TREND</span>
          <h3>Linear in adults</h3>
          <MiniPlot kind="adult" />
          <p>Does expression tend to rise or fall with age among adults? This model summarizes the association with one average rate of change across adult ages, beginning at 20.</p>
          <p className="age-model-takeaway"><strong>Useful for:</strong> finding a steady overall adult-age association.</p>
        </article>
        <article className="age-model-card">
          <span className="age-model-number">02 / EARLIER LIFE</span>
          <h3>Linear in juveniles</h3>
          <MiniPlot kind="juvenile" />
          <p>Does expression change with age before 20? A shared developmental model estimates a separate rate for younger people and for adults. It also allows the groups to differ at the boundary.</p>
          <p className="age-model-takeaway"><strong>Useful for:</strong> comparing early-life and adult-age patterns when both are supported by the data.</p>
        </article>
        <article className="age-model-card">
          <span className="age-model-number">03 / AGE GROUPS</span>
          <h3>Age-bin model</h3>
          <MiniPlot kind="bins" />
          <p>How does each age group compare with adults aged 20–30? Supported ten-year groups are compared with that reference, allowing the pattern to change direction rather than follow one straight line.</p>
          <p className="age-model-takeaway"><strong>Useful for:</strong> seeing patterns that a single average trend can miss.</p>
        </article>
      </div>
      <p className="age-model-note">The adult trend and the adult slope in the developmental model come from different fits. Availability depends on the samples and age groups represented in each analysis. The age ranges shown reflect the default configuration.</p>
    </section>
  )
}
