import { useEffect, useMemo, useState } from 'react'
import { assetFailureKind, type AssetFailureKind } from '../../data/release'
import {
  loadCrossType,
  loadCrossTypeGene,
  type CrossTypeIndex,
  type CrossTypeGene,
} from '../../data/crossType'
import { formatNumber } from '../../theme/format'
import { orderedResults, patternExtent, patternPoints, visibleGenes } from './crossTypeView'
import './crossType.css'

export function CrossTypePage({
  preferredGene,
  onSelectGene,
  onOpenGene,
  onOpenCellType,
  onOpenTrajectories,
}: {
  preferredGene: string | null
  onSelectGene: (gene: string) => void
  onOpenGene: (gene: string) => void
  onOpenCellType: (gene: string, cellType: string) => void
  onOpenTrajectories: () => void
}) {
  const [index, setIndex] = useState<CrossTypeIndex | null>(null)
  const [error, setError] = useState<AssetFailureKind | null>(null)
  const [geneError, setGeneError] = useState(false)
  const [detail, setDetail] = useState<CrossTypeGene | null>(null)
  const [query, setQuery] = useState('')
  const [minimum, setMinimum] = useState(0)
  const [limit, setLimit] = useState(30)
  const [selected, setSelected] = useState<string | null>(null)
  useEffect(() => {
    let active = true
    loadCrossType()
      .then((value) => {
        if (active) setIndex(value)
      })
      .catch((failure) => {
        if (active) setError(assetFailureKind(failure))
      })
    return () => {
      active = false
    }
  }, [])
  const gene =
    index?.genes.find((row) => row.gene === preferredGene)?.gene ?? selected ?? index?.genes[0]?.gene ?? null
  useEffect(() => {
    if (!index || !gene) return
    let active = true
    setDetail(null)
    setGeneError(false)
    loadCrossTypeGene(index, gene)
      .then((value) => {
        if (active) setDetail(value)
      })
      .catch(() => {
        if (active) setGeneError(true)
      })
    return () => {
      active = false
    }
  }, [index, gene])
  const visible = useMemo(() => (index ? visibleGenes(index, query, minimum) : []), [index, query, minimum])
  const choose = (value: string) => {
    setSelected(value)
    onSelectGene(value)
  }
  const percent = (value: number | null) => (value === null ? 'Unavailable' : `${formatNumber(value * 100)}%`)
  if (!index)
    return (
      <section className="content-page cross-type-page">
        <h1>Across cell types</h1>
        <p role="status">
          {error === 'unavailable'
            ? 'Cross-type summaries are not available in this release.'
            : error
              ? 'Cross-type summaries could not be loaded for this release.'
              : 'Loading cross-type summaries…'}
        </p>
      </section>
    )
  const extent = patternExtent(index.clusters)
  return (
    <section className="content-page cross-type-page">
      <header className="cross-type-heading">
        <div>
          <p className="eyebrow">Shared age-associated expression</p>
          <h1>Across cell types</h1>
          <p>Which genes change with age in several cell types, and do their age-bin patterns agree?</p>
        </div>
        <button onClick={onOpenTrajectories}>Explore age-bin trajectories</button>
      </header>
      <p className="cross-type-scope">
        Each cell type has its own omnibus age-bin test. Recurrence counts significant cell types; it does not
        measure replication across studies or imply matching effect directions. Study support remains
        available in the gene and cell-type views.
      </p>
      <div className="cross-type-stats" aria-label="Cross-type analysis coverage">
        <div>
          <strong>{index.cellTypes.length}</strong>
          <span>cell types with results</span>
        </div>
        <div>
          <strong>{index.genes.length.toLocaleString()}</strong>
          <span>genes represented</span>
        </div>
        <div>
          <strong>{index.genes.filter((row) => row.significant > 1).length.toLocaleString()}</strong>
          <span>significant in multiple types</span>
        </div>
      </div>
      <section aria-labelledby="recurrence-heading">
        <h2 id="recurrence-heading">Shared age-bin signals</h2>
        <p>
          Significance: omnibus FDR &lt; {formatNumber(index.deFdr)}. The denominator is the number of cell
          types with a finite test for that gene.
        </p>
        <div className="cross-type-controls">
          <label>
            Find a gene
            <input
              type="search"
              value={query}
              onChange={(event) => {
                setQuery(event.target.value)
                setLimit(30)
              }}
            />
          </label>
          <label>
            Minimum significant cell types
            <select
              value={minimum}
              onChange={(event) => {
                setMinimum(Number(event.target.value))
                setLimit(30)
              }}
            >
              {Array.from({ length: index.cellTypes.length + 1 }, (_, value) => (
                <option key={value} value={value}>
                  {value === 0 ? 'Any' : value}
                </option>
              ))}
            </select>
          </label>
        </div>
        <div className="cross-type-layout">
          <div className="cross-type-gene-list">
            <table>
              <thead>
                <tr>
                  <th>Gene</th>
                  <th>Significant / tested types</th>
                </tr>
              </thead>
              <tbody>
                {visible.slice(0, limit).map((row) => (
                  <tr key={row.gene} className={row.gene === gene ? 'is-selected' : ''}>
                    <th>
                      <button aria-pressed={row.gene === gene} onClick={() => choose(row.gene)}>
                        {row.gene}
                      </button>
                    </th>
                    <td>
                      {row.significant} / {row.tested}
                    </td>
                  </tr>
                ))}
              </tbody>
            </table>
            {!visible.length && <p role="status">No genes match these filters.</p>}
            {visible.length > limit && <button onClick={() => setLimit(limit + 30)}>Show more genes</button>}
          </div>
          <section className="cross-type-detail" aria-label="Selected gene across cell types">
            <h3>{gene ?? 'Select a gene'}</h3>
            {gene && <button onClick={() => onOpenGene(gene)}>Open gene view</button>}
            {geneError ? (
              <p role="status">Results for this gene could not be loaded.</p>
            ) : !detail || detail.gene !== gene ? (
              <p role="status">Loading gene results…</p>
            ) : (
              <>
                <p>
                  Largest age-bin effect relative to {index.referenceBin} years. Positive and negative values
                  refer to different age-bin contrasts; inspect the bin before comparing magnitudes.
                </p>
                <div className="cross-type-table">
                  <table>
                    <thead>
                      <tr>
                        <th>Cell type</th>
                        <th>Age-bin test FDR</th>
                        <th>Largest effect (log₂ FC)</th>
                        <th>Age bin</th>
                      </tr>
                    </thead>
                    <tbody>
                      {orderedResults(index, detail).map(({ cellType, result }) => (
                        <tr key={cellType}>
                          <th>
                            <button onClick={() => onOpenCellType(detail.gene, cellType)}>{cellType}</button>
                          </th>
                          <td>
                            {!result || !result.tested ? (
                              'Not tested'
                            ) : (
                              <span className={result.significant ? 'cross-type-significant' : ''}>
                                {formatNumber(result.padj!)}
                                {result.significant ? ' · significant' : ' · not significant'}
                              </span>
                            )}
                          </td>
                          <td
                            className={
                              result?.effect == null
                                ? ''
                                : result.effect >= 0
                                  ? 'cross-type-positive'
                                  : 'cross-type-negative'
                            }
                          >
                            {result?.effect == null ? 'Unavailable' : formatNumber(result.effect)}
                          </td>
                          <td>{result?.ageBin ?? 'Unavailable'}</td>
                        </tr>
                      ))}
                    </tbody>
                  </table>
                </div>
                <h4>Trajectory agreement</h4>
                {detail.pattern ? (
                  <>
                    <p>
                      Among {detail.pattern.types} cell types with eligible trajectories at FDR ≤{' '}
                      {formatNumber(index.clusterFdr)}:
                    </p>
                    <dl className="cross-type-agreement">
                      <div>
                        <dt>In the most common pattern</dt>
                        <dd>{percent(detail.pattern.dominantFraction)}</dd>
                      </div>
                      <div>
                        <dt>Pairs in the same pattern</dt>
                        <dd>{percent(detail.pattern.sameClusterFraction)}</dd>
                      </div>
                      <div>
                        <dt>Mean pairwise shape distance</dt>
                        <dd>
                          {detail.pattern.meanDistance === null
                            ? 'Unavailable'
                            : formatNumber(detail.pattern.meanDistance)}
                        </dd>
                      </div>
                    </dl>
                    <p className="cross-type-note">
                      Lower distance means more similar standardized age-bin profiles. This describes shape
                      agreement, not independent replication.
                    </p>
                  </>
                ) : (
                  <p>No trajectory-agreement summary is available for this gene.</p>
                )}
              </>
            )}
          </section>
        </div>
      </section>
      <section aria-labelledby="patterns-heading">
        <h2 id="patterns-heading">Shared trajectory patterns</h2>
        <p>
          Mean standardized shapes from the combined trajectory analysis. Scaling emphasizes shape rather than
          effect size. All patterns share the same vertical scale.
        </p>
        {!index.clusters.length ? (
          <p>No shared trajectory patterns meet the support requirements.</p>
        ) : (
          <div className="cross-type-patterns">
            {index.clusters.map((cluster) => (
              <figure key={cluster.id}>
                <figcaption>
                  <strong>Pattern {cluster.id}</strong>
                  <span>{cluster.count} gene–cell-type profiles</span>
                </figcaption>
                <svg
                  viewBox="0 0 320 150"
                  role="img"
                  aria-label={`Mean standardized trajectory for pattern ${cluster.id}`}
                >
                  <line x1="30" x2="290" y1="78" y2="78" stroke="currentColor" opacity=".2" />
                  <text x="25" y="31" textAnchor="end" fontSize="9" data-axis="pattern-extent">
                    {formatNumber(extent)}
                  </text>
                  <text x="25" y="81" textAnchor="end" fontSize="9">
                    0
                  </text>
                  <text x="25" y="131" textAnchor="end" fontSize="9">
                    {formatNumber(-extent)}
                  </text>
                  <polyline
                    points={patternPoints(cluster.mean, extent)}
                    fill="none"
                    stroke="var(--accent, #287d70)"
                    strokeWidth="3"
                  />
                  <text x="30" y="145" fontSize="11">
                    {index.commonBins[0]}
                  </text>
                  <text x="290" y="145" textAnchor="end" fontSize="11">
                    {index.commonBins.at(-1)}
                  </text>
                </svg>
              </figure>
            ))}
          </div>
        )}
      </section>
      <details>
        <summary>Age support and excluded cell types</summary>
        <p>
          Shared bins: {index.commonBins.join(', ') || 'none'} years. Agreement summaries use the stricter
          trajectory threshold and shared-bin support.
        </p>
        {index.exclusions.length ? (
          <ul>
            {index.exclusions.map((row) => (
              <li key={row.cellType}>
                <strong>{row.cellType}:</strong> {row.reason}
              </li>
            ))}
          </ul>
        ) : (
          <p>No exclusions were recorded for the shared-bin analysis.</p>
        )}
      </details>
    </section>
  )
}
