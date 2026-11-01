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
  const { context, analysis, provenance } = release
  return (
    <>
      <h2>Analysis and provenance</h2>
      <dl>
        <dt>Analysis generated</dt>
        <dd>{analysis?.timestamp || 'Not recorded'}</dd>
        <dt>Method</dt>
        <dd>
          {analysis?.method
            ? `${analysis.method}${analysis.methodVersion ? ` ${analysis.methodVersion}` : ''}`
            : 'Not recorded'}
        </dd>
        <dt>Cell-type annotation</dt>
        <dd>{analysis?.splitBy || 'Not recorded'}</dd>
        <dt>Analysis configuration SHA-256</dt>
        <dd>{analysis?.configurationSha256 || 'Not recorded'}</dd>
        <dt>Model metadata comparison</dt>
        <dd>
          {provenance?.status || 'Not checked'}.{' '}
          {provenance?.note || 'This release does not contain a metadata comparison.'}
        </dd>
      </dl>
      <p>
        The export date and ID identify this site export. They do not identify the analysis run. Source files
        are fingerprinted in the manifest, but existing reports and tables do not record enough lineage to
        prove they came from the same analysis.
      </p>
      <h2>Studies in this release</h2>
      {context?.studies.length ? (
        <>
          <p>
            Study identifiers come from this release. Citation and annotation details are from the study
            configuration.
          </p>
          <ul>
            {context.studies.map((study) => (
              <li key={study.id}>
                <h3>{study.id}</h3>
                {study.inModels !== undefined && (
                  <p>
                    {study.harmonizationAvailable
                      ? 'Harmonization diagnostics available.'
                      : 'No harmonization diagnostics in this release.'}{' '}
                    {study.inModels
                      ? 'Represented in exported DE models.'
                      : 'Not represented in exported DE models.'}
                  </p>
                )}
                <p>
                  <Doi value={study.publicationDoi} />
                </p>
                {study.annotation && <p>Cell types: {study.annotation}</p>}
                {study.technology && <p>Technology: {study.technology}</p>}
                {study.inclusion && <p>Model inclusion: {study.inclusion}</p>}
                {study.siteComments && <p>{study.siteComments}</p>}
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
              </li>
            ))}
          </ul>
        </>
      ) : (
        <p>No study identifiers were recorded in the exported models.</p>
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
      <p>
        {release.analysis?.method
          ? `DE method: ${release.analysis.method}${release.analysis.methodVersion ? ` ${release.analysis.methodVersion}` : ''}.`
          : 'DE method metadata is unavailable.'}
      </p>
      {inclusion && (
        <p>
          Configured DE inclusion:{' '}
          {inclusion.minimum_age_years_inclusive != null
            ? `age ≥${inclusion.minimum_age_years_inclusive}`
            : 'minimum age not recorded'}
          ;{' '}
          {inclusion.minimum_cells_filter_enabled === false
            ? 'minimum-cell filter disabled'
            : inclusion.minimum_cells_per_pseudobulk_inclusive != null
              ? `minimum ${inclusion.minimum_cells_per_pseudobulk_inclusive} cells per pseudobulk`
              : 'minimum-cell setting not recorded'}
          . Model-specific sample counts and designs are available in differential expression.
        </p>
      )}
    </>
  )
}

export function SourcePermissions({ release }: { release: ReleaseManifest }) {
  return (
    <>
      <h2>Source and report permissions</h2>
      <p>
        Study licences are not recorded in this release. Review the terms and attribution requirements for
        each source before distribution.
      </p>
      {release.pathwaySources?.length ? (
        <p>
          Pathway source files: {release.pathwaySources.map((source) => source.path).join(', ')}. Their
          fingerprints are in the manifest; licence and attribution review is still required.
        </p>
      ) : (
        <p>No pathway source files were recorded.</p>
      )}
      {release.pathwaySources?.length ? (
        <p>
          For sources from MSigDB, review the{' '}
          <a
            href="https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp"
            target="_blank"
            rel="noreferrer"
          >
            MSigDB licence and source-specific conditions
          </a>
          .
        </p>
      ) : null}
      <p>
        Included HTML reports may contain sample-level identifiers and diagnostic plots. Their companion
        downloads are not included in this export. Review report contents and sharing permissions for the
        intended audience.
      </p>
    </>
  )
}
