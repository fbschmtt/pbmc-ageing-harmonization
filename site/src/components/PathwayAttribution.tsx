import { useId, useState } from 'react'
import type { PathwaySet } from '../data/pathways'

function isWebLink(value: string) {
  try {
    const url = new URL(value)
    return url.protocol === 'https:' || url.protocol === 'http:'
  } catch {
    return false
  }
}

export function PathwayAttribution({ activePathway }: { activePathway: PathwaySet | null }) {
  const [isOpen, setIsOpen] = useState(false)
  const panelId = useId()

  return (
    <div
      className="pathway-attribution"
      onBlur={(event) => {
        if (!event.currentTarget.contains(event.relatedTarget)) setIsOpen(false)
      }}
      onKeyDown={(event) => {
        if (event.key === 'Escape') setIsOpen(false)
      }}
    >
      <button
        type="button"
        className="pathway-attribution-trigger"
        title="Show pathway sources, licence, and attribution"
        aria-expanded={isOpen}
        aria-controls={panelId}
        onClick={() => setIsOpen((open) => !open)}
      >
        Licence &amp; attribution
        <span aria-hidden="true">ⓘ</span>
      </button>
      {isOpen && (
        <div id={panelId} className="pathway-attribution-panel">
          <strong>MSigDB gene sets</strong>
          <p>
            The curated sets shown here were obtained from the Molecular Signatures Database. Current MSigDB
            releases are provided under CC BY 4.0, with additional terms for some source collections.
          </p>
          <p>
            Attribute the Broad Institute, MIT, and Regents of the University of California, and review the{' '}
            <a
              href="https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp"
              target="_blank"
              rel="noreferrer"
            >
              complete MSigDB licence terms
            </a>
            .
          </p>
          {activePathway && (
            <p>
              Selected set:{' '}
              {isWebLink(activePathway.description) ? (
                <a href={activePathway.description} target="_blank" rel="noreferrer">
                  {activePathway.id}
                </a>
              ) : (
                activePathway.id
              )}
              <span className="pathway-attribution-file">Source file: {activePathway.sourceFile}</span>
            </p>
          )}
        </div>
      )}
    </div>
  )
}
