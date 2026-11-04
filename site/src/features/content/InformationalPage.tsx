import { useEffect, useState } from 'react'
import { loadReleaseManifest, type ReleaseManifest } from '../../data/release'
import { DatasetContext, ModellingContext, SourcePermissions } from './ReleaseContext'
import type { Page } from '../../app/navigation'
import { AgeModelGuide } from './AgeModelGuide'

const copy: Record<Extract<Page, 'datasets' | 'modelling'>, { title: string; body: string }> = {
  datasets: {
    title: 'Datasets',
    body: 'The atlas brings together human peripheral blood studies with different cohorts, technologies and available metadata.',
  },
  modelling: {
    title: 'Models and interpretation',
    body: 'The analyses compare gene expression and cell-type fractions across samples, with attention to differences between studies and the information available for each model.',
  },
}

export function InformationalPage({ page }: { page: Extract<Page, 'datasets' | 'modelling' | 'download'> }) {
  const [release, setRelease] = useState<ReleaseManifest | null>(null)
  useEffect(() => {
    if (page === 'download') return
    loadReleaseManifest()
      .then(setRelease)
      .catch(() => setRelease(null))
  }, [page])
  const reports = release?.assets.filter((asset) => asset.kind === 'report') ?? []
  const reportIndex = release?.assets.find((asset) => asset.kind === 'report_index')
  const warnings = [...new Set(reports.map((report) => report.warning).filter(Boolean))]
  if (page === 'download')
    return (
      <section className="content-page">
        <p className="eyebrow">Download</p>
        <h1>Downloads</h1>
        <p>Downloads will be available in a future release.</p>
      </section>
    )
  const content = copy[page]
  return (
    <section className="content-page information-page">
      <p className="eyebrow">{page}</p>
      <h1>{content.title}</h1>
      <p>{content.body}</p>
      {page === 'datasets' && (
        <p>
          Cell types are harmonized using the Allen Institute for Immunology (AIFI) reference. The studies
          represented in a model depend on sample and metadata availability; see the model details for its
          coverage.
        </p>
      )}
      {page === 'modelling' && (
        <>
          <h2>Gene expression</h2>
          <p>
            Counts are summed for each sample and cell type, then analysed with per-study and combined
            pseudobulk models. The combined age model adjusts for study site where estimable, sex and
            pseudobulk depth; models for other covariates use samples with the required information.
          </p>
          <p>
            Age-bin models compare supported age groups with ages 20–30. The trajectory views show only genes
            selected for clustering, not all tested genes. Associations do not establish causation.
          </p>
          <AgeModelGuide />
          <h2>Cell composition</h2>
          <p>
            Fraction models are fitted separately for each study. They describe relative representation within
            a parent cell compartment or among retained PBMCs, not absolute cell abundance.
          </p>
        </>
      )}
      {page === 'datasets' && (
        <p>
          <a href={`${import.meta.env.BASE_URL}data/release.json`}>View the release manifest</a>
        </p>
      )}
      {page === 'datasets' && release && <DatasetContext release={release} />}
      {page === 'modelling' && (
        <section className="information-section">
          <h2>Detailed model diagnostics</h2>
          {reportIndex ? (
            <>
              <p>
                Browse the included differential-expression and cross-cell-type reports in the report index.
                Some companion downloads linked inside reports are not included in this site.
              </p>
              {warnings.map((warning) => (
                <p key={warning} role="status">
                  {warning}
                </p>
              ))}
              <p>
                <a
                  href={`${import.meta.env.BASE_URL}data/${reportIndex.path}`}
                  target="_blank"
                  rel="noopener noreferrer"
                >
                  Open the report index
                </a>
              </p>
            </>
          ) : reports.length ? (
            <p>The report index is unavailable for this release.</p>
          ) : (
            <p>No detailed reports are included in this release.</p>
          )}
        </section>
      )}
      {page === 'modelling' && (
        <section className="information-metadata" aria-label="Release details">
          <h2>Release details</h2>
          {release && <ModellingContext release={release} />}
          <p>
            <a href={`${import.meta.env.BASE_URL}data/release.json`}>View the release manifest</a>
          </p>
          {release && <SourcePermissions release={release} />}
        </section>
      )}
    </section>
  )
}
