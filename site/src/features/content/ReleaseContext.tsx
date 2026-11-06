import type { ReleaseManifest } from '../../data/release'

function Doi({ value }: { value: string | null }) {
  if (!value || !/^10\.\d{4,9}\/[^\s]+$/.test(value)) return <>Citation unavailable</>
  return (
    <a href={`https://doi.org/${encodeURI(value)}`} target="_blank" rel="noreferrer">
      {value}
    </a>
  )
}

export function DatasetContext({ release }: { release: ReleaseManifest }) {
  const { context } = release
  return (
    <>
      {context?.studies.length ? (
        <>
          <section className="study-index" aria-labelledby="study-index-heading">
            <h2 id="study-index-heading">Studies in this release</h2>
            <ul>
              {context.studies.map((study, index) => (
                <li key={study.id}>
                  <a
                    href={`#study-section-${index}`}
                    onClick={(event) => {
                      event.preventDefault()
                      document.getElementById(`study-section-${index}`)?.scrollIntoView({
                        behavior: 'auto',
                        block: 'start',
                      })
                    }}
                  >
                    {study.id}
                  </a>
                </li>
              ))}
            </ul>
          </section>
          <div className="dataset-studies">
            {context.studies.map((study, index) => (
              <section id={`study-section-${index}`} key={study.id} className="dataset-study">
                <p className="eyebrow">Study</p>
                <h3>{study.publicationTitle || study.id}</h3>
                {study.publicationCitation && <p className="study-citation">{study.publicationCitation}</p>}
                {study.publicationDoi && <p className="study-doi"><Doi value={study.publicationDoi} /></p>}
                {!study.publicationTitle && !study.publicationDoi && (
                  <p>Publication details have not been provided for this study.</p>
                )}
                <p className="study-counts">
                  {study.donors != null || study.cells != null ? (
                    <>
                      In this release: {study.donors != null && `${study.donors.toLocaleString()} donors`}
                      {study.donors != null && study.cells != null && ' · '}
                      {study.cells != null && `${study.cells.toLocaleString()} cells`}
                    </>
                  ) : (
                    'Donor and cell counts are not available for this study in this release.'
                  )}
                </p>
                <StudyFigures release={release} study={study.id} />
                {release.assets
                  .filter((asset) => asset.kind === 'report' && asset.study === study.id)
                  .map((asset) => (
                    <p key={asset.path}>
                      <a
                        href={`${import.meta.env.BASE_URL}data/${asset.path}`}
                        target="_blank"
                        rel="noreferrer"
                      >
                        Open full harmonization report
                      </a>
                    </p>
                  ))}
              </section>
            ))}
          </div>
        </>
      ) : (
        <p>No study identifiers were recorded for this release.</p>
      )}
    </>
  )
}

function StudyFigures({ release, study }: { release: ReleaseManifest; study: string }) {
  const figures = release.assets.filter(
    (asset) => asset.study === study && (asset.kind === 'study_umap' || asset.kind === 'study_heatmap'),
  )
  return figures.length ? (
    <div className="study-figure-list">
      {(['study_umap', 'study_heatmap'] as const).map((kind) => {
        const figure = figures.find((asset) => asset.kind === kind)
        return figure ? (
          <figure key={kind}>
            <img
              src={`${import.meta.env.BASE_URL}data/${figure.path}`}
              loading="lazy"
              alt={`${study} ${kind === 'study_umap' ? 'UMAP colored by AIFI L2 labels' : 'study cell type versus AIFI L2 heatmap'}`}
            />
            <figcaption>
              {kind === 'study_umap' ? 'Harmonized UMAP · AIFI L2 labels' : 'Study types × AIFI L2'}
            </figcaption>
          </figure>
        ) : null
      })}
    </div>
  ) : null
}

export function ModellingContext({ release }: { release: ReleaseManifest }) {
  const inclusion = release.analysis?.sampleInclusion
  return (
    <>
      {release.analysis?.method && (
        <p>
          DE method: {release.analysis.method}
          {release.analysis.methodVersion ? ` ${release.analysis.methodVersion}` : ''}.
        </p>
      )}
      {inclusion && (
        <p>
          {inclusion.minimum_age_years_inclusive != null &&
            `Samples younger than ${inclusion.minimum_age_years_inclusive} years are excluded. `}
          {inclusion.minimum_cells_filter_enabled === false
            ? 'No minimum cell count is applied to pseudobulk samples.'
            : inclusion.minimum_cells_per_pseudobulk_inclusive != null
              ? `Pseudobulk samples require at least ${inclusion.minimum_cells_per_pseudobulk_inclusive} cells.`
              : null}
        </p>
      )}
    </>
  )
}
