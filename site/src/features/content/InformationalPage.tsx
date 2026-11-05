import { useEffect, useState } from 'react'
import { loadReleaseManifest, type ReleaseManifest } from '../../data/release'
import { DatasetContext, ModellingContext } from './ReleaseContext'
import type { Page } from '../../app/navigation'
import { AgeModelGuide } from './AgeModelGuide'

const copy: Record<Extract<Page, 'datasets' | 'modelling' | 'about'>, { title: string; body: string }> = {
  datasets: {
    title: 'Datasets',
    body: 'Explore the human blood studies behind the atlas and the research teams whose work makes these comparisons possible.',
  },
  modelling: {
    title: 'Models and interpretation',
    body: 'The analyses compare gene expression and cell-type fractions across samples, with attention to differences between studies and the information available for each model.',
  },
  about: {
    title: 'About Blood Ageing Explorer',
    body: 'Blood Ageing Explorer brings together published single-cell studies of human peripheral blood to make age-related patterns across immune cell types easier to explore.',
  },
}

export function InformationalPage({
  page,
  onNavigate,
}: {
  page: Extract<Page, 'datasets' | 'modelling' | 'about'>
  onNavigate?: (page: Page) => void
}) {
  const [release, setRelease] = useState<ReleaseManifest | null>(null)
  useEffect(() => {
    if (page === 'about') return
    loadReleaseManifest()
      .then(setRelease)
      .catch(() => setRelease(null))
  }, [page])
  const reports = release?.assets.filter((asset) => asset.kind === 'report') ?? []
  const reportIndex = release?.assets.find((asset) => asset.kind === 'report_index')
  const warnings = [...new Set(reports.map((report) => report.warning).filter(Boolean))]
  const content = copy[page]
  return (
    <section className="content-page information-page">
      <p className="eyebrow">{page}</p>
      <h1>{content.title}</h1>
      <p>{content.body}</p>
      {page === 'datasets' && (
        <p>
          This atlas builds on the work of the teams who collected and shared these cohorts. Browse the
          publications below to see each study’s contribution and source paper.
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
          <h2>Cell states and age-associated expression</h2>
          <p>
            One question is whether age-associated differential expression mostly reflects shifts in cell states
            or more independent changes in individual genes—in other words, how coordinated those changes are.
            The <a href="#singleCell" onClick={(event) => {
              event.preventDefault()
              onNavigate?.('singleCell')
            }}>Cell states view</a> lets you explore local clusters and age-associated genes within each cell type.
          </p>
          <AgeModelGuide />
          <h2>Cell composition</h2>
          <p>
            Fraction models are fitted separately for each study. They describe relative representation within
            a parent cell compartment or among retained PBMCs, not absolute cell abundance.
          </p>
        </>
      )}
      {page === 'datasets' && release && <DatasetContext release={release} />}
      {page === 'about' && (
        <section className="information-section">
          <h2>Purpose</h2>
          <p>
            The explorer presents summary analyses of published human peripheral blood datasets. It is
            intended for browsing and generating research questions; observed associations are not causal
            conclusions or clinical guidance.
          </p>
          <h2>Studies and methods</h2>
          <p>
            Study sources and publication records are listed on the Datasets page. The Modelling page
            describes how expression and cell-composition results are estimated and how to interpret them.
          </p>
        </section>
      )}
      {page === 'modelling' && (
        <section className="information-section">
          <h2>Detailed model diagnostics</h2>
          {reportIndex ? (
            <>
              <p>
                Browse the included differential-expression and cross-cell-type reports in the report index.
                The reports include detailed model results and diagnostic plots.
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
          ) : release?.assets.some((asset) => asset.kind === 'report') ? (
            <p>The report index is unavailable for this release.</p>
          ) : (
            <p>No detailed reports are included in this release.</p>
          )}
        </section>
      )}
      {page === 'modelling' && (
        <section className="information-metadata" aria-label="Release details">
          <h2>Analysis details</h2>
          {release && <ModellingContext release={release} />}
        </section>
      )}
    </section>
  )
}
