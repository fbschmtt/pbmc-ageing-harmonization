import type { TrajectoryRecord } from '../data/atlas'
import { formatCpm, formatNumber } from '../theme/format'
import { TrajectoryChart } from './TrajectoryChart'
import { useEffect, useState } from 'react'
import { loadCombinedAgeEffect } from '../data/differential'
import './geneTypePairCard.css'

export type GeneTypePairActions = {
  onOpenGene?: (record: TrajectoryRecord) => void
  onOpenCellType?: (record: TrajectoryRecord) => void
  onOpenTrajectory?: (record: TrajectoryRecord) => void
}

export function GeneTypePairCard({
  record,
  actions,
}: {
  record: TrajectoryRecord
  actions: GeneTypePairActions
}) {
  const [combinedAgeEffect, setCombinedAgeEffect] = useState<number | null>(null)
  const [ageEffectLoaded, setAgeEffectLoaded] = useState(false)
  useEffect(() => {
    let active = true
    setAgeEffectLoaded(false)
    loadCombinedAgeEffect(record.cellType, record.gene)
      .then((value) => {
        if (active) setCombinedAgeEffect(value)
      })
      .catch(() => {
        if (active) setCombinedAgeEffect(null)
      })
      .finally(() => {
        if (active) setAgeEffectLoaded(true)
      })
    return () => { active = false }
  }, [record.cellType, record.gene])
  return (
    <article className="pair-card">
      <p className="eyebrow">Gene–celltype pair</p>
      <h2>
        {actions.onOpenGene ? (
          <button className="pair-card-title-link" onClick={() => actions.onOpenGene?.(record)}>
            {record.gene}
          </button>
        ) : (
          record.gene
        )}
      </h2>
      <p className="detail-subtitle">
        {actions.onOpenCellType ? (
          <button className="pair-card-cell-link" onClick={() => actions.onOpenCellType?.(record)}>
            {record.cellType}
          </button>
        ) : (
          record.cellType
        )}{' '}
        · {record.tissue}
      </p>
      <div className="trajectory-figure">
        <div className="trajectory-chart-title">Scaled age-bin effect</div>
        <TrajectoryChart record={record} />
        <div className="trajectory-labels">
          <span>20–30</span>
          <span>50–60</span>
          <span>80–90</span>
        </div>
      </div>
      <dl>
        <div>
          <dt>Reference expression (20–30)</dt>
          <dd>
            {record.referenceExpression != null ? formatCpm(record.referenceExpression) : 'Unavailable'}
          </dd>
        </div>
        <div>
          <dt>Age-bin omnibus adjusted p-value</dt>
          <dd>{formatNumber(record.adjustedPValue)}</dd>
        </div>
        <div>
          <dt>Combined linear age effect per decade</dt>
          <dd>
            {ageEffectLoaded && combinedAgeEffect !== null
              ? `${combinedAgeEffect * 10 > 0 ? '+' : ''}${formatNumber(combinedAgeEffect * 10)} log₂ fold change`
              : ageEffectLoaded ? 'Unavailable' : 'Loading…'}
          </dd>
        </div>
        <div>
          <dt>Signed maximum age-bin change (vs 20–30)</dt>
          <dd>
            {record.valueAtMaxAbsChange > 0 ? '+' : ''}{formatNumber(record.valueAtMaxAbsChange)} log₂ fold change
          </dd>
        </div>
        <div>
          <dt>Shape cluster</dt>
          <dd>{record.cluster}</dd>
        </div>
        <div>
          <dt>Mean cells per donor</dt>
          <dd>{formatNumber(record.meanCellsPerDonor)}</dd>
        </div>
      </dl>
      {(actions.onOpenGene || actions.onOpenCellType || actions.onOpenTrajectory) && (
        <div className="pair-card-actions">
          {actions.onOpenGene && (
            <button className="text-button" onClick={() => actions.onOpenGene?.(record)}>
              View gene profile →
            </button>
          )}
          {actions.onOpenCellType && (
            <button className="text-button" onClick={() => actions.onOpenCellType?.(record)}>
              View differential expression →
            </button>
          )}
          {actions.onOpenTrajectory && (
            <button className="text-button" onClick={() => actions.onOpenTrajectory?.(record)}>
              Open age-bin trajectory map →
            </button>
          )}
        </div>
      )}
    </article>
  )
}
