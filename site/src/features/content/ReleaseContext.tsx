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
      <h2>Studies represented in exported models</h2>
      {context?.studies.length ? (
        <>
          <p>
            Study identifiers come from the exported models. Citation and annotation details below are from
            the exporter checkout's study configuration; they have not been verified against an analysis-time
            configuration snapshot.
          </p>
          <ul>
            {context.studies.map((study) => (
              <li key={study.id}>
                <strong>{study.id}</strong> · <Doi value={study.publicationDoi} />
                {study.annotation && <p>Cell types: {study.annotation}</p>}
                {study.technology && <p>Technology: {study.technology}</p>}
                {study.inclusion && <p>Model inclusion: {study.inclusion}</p>}
                {study.siteComments && <p>{study.siteComments}</p>}
                {release.assets
                  .filter((asset) => asset.kind === 'report' && asset.study === study.id)
                  .map((asset) => (
                    <p key={asset.path}>
                      <a
                        href={`${import.meta.env.BASE_URL}data/${asset.path}`}
                        target="_blank"
                        rel="noreferrer"
                      >
                        Open harmonization diagnostics
                      </a>
                    </p>
                  ))}
                <p>
                  Data licence: {study.licence || 'Not recorded; review source terms before publication.'}
                </p>
              </li>
            ))}
          </ul>
        </>
      ) : (
        <p>No study identifiers were recorded in the exported models.</p>
      )}
      <h2>Age-bin support</h2>
      {context?.ageSupport?.commonBins?.length ? (
        <>
          <p>
            Common cross-type trajectory bins: {context.ageSupport.commonBins.join(', ')}.{' '}
            {context.ageSupport.minimumSamplesPerBin != null &&
              `Configured minimum: ${context.ageSupport.minimumSamplesPerBin} samples per bin.`}
          </p>
          {context.ageSupport.cellTypes.length > 0 && (
            <details>
              <summary>Cell-type support and exclusions (sample counts per bin)</summary>
              <ul>
                {context.ageSupport.cellTypes.map((cell) => (
                  <li key={cell.cellType}>
                    <strong>{cell.cellType}</strong>: {cell.status}
                    {cell.reason && ` — ${cell.reason}`}
                    {cell.ageBinCounts && (
                      <span>
                        {' '}
                        ·{' '}
                        {Object.entries(cell.ageBinCounts)
                          .map(([bin, count]) => `${bin}: ${count}`)
                          .join(', ')}
                      </span>
                    )}
                  </li>
                ))}
              </ul>
            </details>
          )}
        </>
      ) : (
        <p>Age-bin support metadata is unavailable for this export.</p>
      )}
    </>
  )
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
