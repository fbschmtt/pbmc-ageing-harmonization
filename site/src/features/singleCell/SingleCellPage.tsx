import { useEffect, useState } from 'react'
import { CellTypeSelector } from '../../components/CellTypeSelector'
import { ChoiceMenu } from '../../components/ChoiceMenu'
import { validAsset, type SingleCell, type SingleCellIndex } from '../../data/contracts.mjs'
import { fetchReleaseAsset } from '../../data/release'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import { formatNumber } from '../../theme/format'
import './singleCell.css'

export function SingleCellPage({
  preferredCellType,
  onSelectCellType,
}: {
  preferredCellType: string | null
  onSelectCellType: (value: string) => void
}) {
  const [index, setIndex] = useState<SingleCellIndex | null>(null)
  const [data, setData] = useState<SingleCell | null>(null)
  const [error, setError] = useState(false)
  const [loading, setLoading] = useState(true)
  const [cluster, setCluster] = useState<string | null>(null)

  const hasPreferred = index?.cellTypes.some((item) => item.label === preferredCellType) ?? false
  const selected = index?.cellTypes.find((item) => item.label === preferredCellType) ?? index?.cellTypes[0]
  const selectedCluster = data?.groups.find((group) => group.cluster === cluster) ?? data?.groups[0]
  const fractionDetail = data?.fractionDetails?.find((item) => item.cluster === selectedCluster?.cluster)
  const imageSource = (image: string) => image.startsWith('single-cell-images/')
    ? `${import.meta.env.BASE_URL}data/${image}?release=${encodeURIComponent(data?.releaseId || '')}`
    : `data:image/png;base64,${image}`
  const maxFitRows = Math.max(0, ...(data?.fractionDetails ?? []).map((item) => item.fits.length))

  useEffect(() => {
    let active = true
    fetchReleaseAsset('single-cell-index.json')
      .then((value) => {
        if (!validAsset('singleCellIndex', value)) throw new Error('Invalid single-cell index')
        if (active) {
          setIndex(value)
          setLoading(false)
        }
      })
      .catch(() => {
        if (active) {
          setError(true)
          setLoading(false)
        }
      })
    return () => {
      active = false
    }
  }, [])

  useEffect(() => {
    if (!selected) return
    let active = true
    setData(null)
    setCluster(null)
    setError(false)
    fetchReleaseAsset(selected.path)
      .then((value) => {
        if (!validAsset('singleCell', value) || value.cellType !== selected.label)
          throw new Error('Invalid single-cell data')
        if (active) setData(value)
      })
      .catch(() => {
        if (active) setError(true)
      })
    return () => {
      active = false
    }
  }, [selected?.path, selected?.label])

  return (
    <section className="single-cell-page">
      <h1 className="visually-hidden">Cell states</h1>
      <header className="single-cell-heading">
        <p className="single-cell-intro">
          Explore local expression patterns, cluster markers, and sample-level cluster fractions.
        </p>
        {index && selected && (
          <div className="single-cell-controls">
            <CellTypeSelector
              label="Cell type"
              value={selected.label}
              groups={groupCellTypes(index.cellTypes.map((item) => item.label)).map((group) => ({
                label: group.label,
                options: group.cellTypes.map((type) => ({
                  value: type,
                  label: type,
                })),
              }))}
              onChange={onSelectCellType}
            />
            {data && (
              <div className="single-cell-stats" aria-label="Selected cell type summary">
                <div>
                  <strong>{data.nCells.toLocaleString()}</strong>
                  <span>cells</span>
                </div>
                <div>
                  <strong>{data.nClusters}</strong>
                  <span>clusters</span>
                </div>
                <div>
                  <strong>{data.integration === 'none' ? 'None' : data.integration}</strong>
                  <span>integration</span>
                </div>
              </div>
            )}
          </div>
        )}
      </header>

      {loading && <p role="status">Loading single-cell analyses…</p>}
      {!loading && !index && !error && <p>Single-cell summaries are unavailable in this release.</p>}
      {index && selected && preferredCellType && !hasPreferred && (
        <p className="single-cell-notice" role="status">
          Cell-state results are unavailable for {preferredCellType}. Showing {selected.label}, the first
          cell type with a completed single-cell analysis.{' '}
          <button type="button" onClick={() => onSelectCellType(selected.label)}>
            Use {selected.label} across views
          </button>
        </p>
      )}
      {error && (
        <p role="status">
          {index
            ? 'The selected single-cell summary could not be loaded.'
            : 'Single-cell summaries could not be loaded.'}
        </p>
      )}
      {index && !index.cellTypes.length && (
        <p>There are no completed single-cell analyses in this release.</p>
      )}
      {index && selected && !data && !error && (
        <p role="status">Loading {selected.label}…</p>
      )}

      {data && (
        <div className="single-cell-layout">
          <article className="single-cell-panel atlas-card">
            <div className="single-cell-panel-heading">
              <div>
                <p className="eyebrow">Local embedding</p>
                <h2>Clusters in context</h2>
              </div>
              <span className="single-cell-panel-meta">{data.nCells.toLocaleString()} cells</span>
            </div>
            <img
              className="single-cell-umap"
              src={imageSource(data.image)}
              alt={`UMAP of ${data.cellType}, colored by local cluster`}
            />
          </article>

          <article className="single-cell-panel atlas-card">
            <div className="single-cell-panel-heading">
              <div>
                <p className="eyebrow">Cluster markers</p>
                <h2>Genes distinguishing clusters</h2>
              </div>
            </div>
            <p>
              Genes with the largest positive and negative log₂ fold changes compared with other clusters in
              this cell type.
            </p>
            {selectedCluster ? (
              <>
                <ChoiceMenu
                  label="Cluster"
                  value={selectedCluster.cluster}
                  groups={[
                    {
                      options: data.groups.map((group) => ({
                        value: group.cluster,
                        label: `Cluster ${group.cluster}`,
                      })),
                    },
                  ]}
                  onChange={setCluster}
                />
                {selectedCluster.higher && selectedCluster.lower ? (
                  <div className="single-cell-clusters">
                    {(
                      [
                        ['Higher in cluster', selectedCluster.higher],
                        ['Lower in cluster', selectedCluster.lower],
                      ] as const
                    ).map(([label, markers]) => (
                      <section key={label} className="single-cell-cluster">
                        <h3>{label}</h3>
                        <div className="single-cell-marker-units">
                          Gene <span>log₂ FC</span>
                        </div>
                        {markers.length ? (
                          <ol>
                            {markers.map((marker) => (
                              <li key={marker.gene}>
                                <strong title={marker.gene}>{marker.gene}</strong>
                                <span>
                                  {marker.effect > 0 ? '+' : ''}
                                  {formatNumber(marker.effect)}
                                </span>
                              </li>
                            ))}
                          </ol>
                        ) : (
                          <p>No markers available.</p>
                        )}
                      </section>
                    ))}
                  </div>
                ) : (
                  <p className="single-cell-empty">Directional markers are unavailable for this cell type.</p>
                )}
              </>
            ) : (
              <p className="single-cell-empty">No cluster markers are available for this analysis.</p>
            )}
            <p className="single-cell-caption">
              Markers describe differences within this cell type, not age associations or causal effects.
            </p>
          </article>

          <article className="single-cell-panel atlas-card">
            <div className="single-cell-panel-heading">
              <div>
                <p className="eyebrow">Age-associated expression</p>
                <h2>Age-DE genes on UMAP</h2>
              </div>
            </div>
            {data.ageDeImage ? (
              <img
                className="single-cell-umap single-cell-age-de-image"
                src={imageSource(data.ageDeImage)}
                alt={`Expression of up to nine significant age-associated genes on the ${data.cellType} UMAP`}
              />
            ) : (
              <p className="single-cell-empty">
                {data.ageDeNote ??
                  'No age-associated expression map was exported for this cell type. It requires merged age-DE results with significant genes detected in this cell type.'}
              </p>
            )}
            <p className="single-cell-caption">
              Expression is shown per cell for genes selected by the sample-level pseudobulk age model; it
              does not imply cell-level significance.
            </p>
          </article>

          <article className="single-cell-panel atlas-card single-cell-fractions">
            <div className="single-cell-panel-heading">
              <div>
                <p className="eyebrow">Cluster composition</p>
                <h2>
                  {selectedCluster
                    ? `Cluster ${selectedCluster.cluster}: fraction across age`
                    : 'Cluster fractions across age'}
                </h2>
              </div>
            </div>
            <p>
              Each point is a study and sample. Fractions are within the selected cell type.
              {data.fractionHasStudyTrends
                ? ' Dashed lines show study trends; the solid line shows the pooled, unadjusted trend.'
                : ' The line shows the pooled, unadjusted trend.'}{' '}
              Lines require samples at different ages.{' '}
              {data.fractionIncludesZeros
                ? 'Eligible samples with no cells in a cluster are included at zero.'
                : 'Samples with no cells in a cluster are omitted.'}
            </p>
            {fractionDetail ? (
              <>
                <div className="single-cell-plot-frame">
                  <img
                    src={imageSource(fractionDetail.image)}
                    alt={`Cluster ${fractionDetail.cluster} sample-level fractions across age for ${data.cellType}, colored by study`}
                  />
                </div>
                <div className="single-cell-fit-scroll">
                  <table className="single-cell-fit-table">
                    <caption>Observed fractions and unadjusted linear fits</caption>
                    <thead>
                      <tr>
                        <th scope="col">Group</th>
                        <th scope="col">Samples</th>
                        <th scope="col">Mean fraction</th>
                        <th scope="col">Slope (pp / decade)</th>
                      </tr>
                    </thead>
                    <tbody>
                      {fractionDetail.fits.map((fit) => (
                        <tr key={fit.study ?? '__pooled__'}>
                          <th scope="row">{fit.study ?? 'Pooled'}</th>
                          <td>{fit.samples}</td>
                          <td>{(fit.meanFraction * 100).toFixed(1)}%</td>
                          <td>
                            {fit.slopePpPerDecade == null
                              ? '—'
                              : `${fit.slopePpPerDecade > 0 ? '+' : ''}${formatNumber(fit.slopePpPerDecade)}`}
                          </td>
                        </tr>
                      ))}
                      {Array.from(
                        { length: maxFitRows - fractionDetail.fits.length },
                        (_, index) => (
                          <tr key={`empty-${index}`} className="single-cell-fit-placeholder" aria-hidden="true">
                            <td>&nbsp;</td>
                            <td>&nbsp;</td>
                            <td>&nbsp;</td>
                            <td>&nbsp;</td>
                          </tr>
                        ),
                      )}
                    </tbody>
                  </table>
                </div>
              </>
            ) : (
              <p className="single-cell-empty">
                {data.fractionNote ??
                  "No cluster-fraction plot was exported. Samples must meet the report's minimum cell-count requirement."}
              </p>
            )}
            {data.fractionNote && fractionDetail && (
              <p className="single-cell-caption">{data.fractionNote}</p>
            )}
          </article>
        </div>
      )}
    </section>
  )
}
