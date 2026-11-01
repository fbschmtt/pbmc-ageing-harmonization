import type { GeneMetadata } from '../data/geneMetadata'

const sections = [
  ['Biological process', 'biologicalProcesses'],
  ['Molecular function', 'molecularFunctions'],
  ['Cellular location', 'cellularComponents'],
] as const

export function GeneFunctionSummary({ metadata, release }: { metadata?: GeneMetadata; release?: string }) {
  return (
    <section className="gene-function-summary" aria-labelledby="gene-function-title">
      <p className="eyebrow" id="gene-function-title">
        Gene function
      </p>
      <dl>
        {sections.map(([label, field]) => {
          const terms = metadata?.goTerms[field] ?? []
          return (
            <div key={field}>
              <dt>{label}</dt>
              <dd title={terms.map((term) => `${term.label} (${term.id})`).join('; ')}>
                {terms.length ? terms.map((term) => term.label).join('; ') : 'No annotation available'}
              </dd>
            </div>
          )
        })}
      </dl>
      <small>Selected GO annotations · {release ?? 'metadata unavailable'}</small>
    </section>
  )
}
