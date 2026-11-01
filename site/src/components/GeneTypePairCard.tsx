import type { TrajectoryRecord } from '../data/atlas'
import type { GeneMetadata } from '../data/geneMetadata'
import { TrajectoryChart } from './TrajectoryChart'
import './geneTypePairCard.css'

export type GeneTypePairActions = {
  onOpenGene?: (record: TrajectoryRecord) => void
  onOpenCellType?: (record: TrajectoryRecord) => void
  onOpenTrajectory?: (record: TrajectoryRecord) => void
}

function formatP(value: number) {
  return value < 0.001 ? value.toExponential(1) : value.toFixed(3)
}

export function GeneTypePairCard({
  record,
  actions,
  metadata,
  metadataRelease,
}: {
  record: TrajectoryRecord
  actions: GeneTypePairActions
  metadata?: GeneMetadata
  metadataRelease?: string
}) {
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
          <dt>Reference expression (20–30, study-balanced CPM)</dt>
          <dd>
            {record.referenceExpression != null ? `${record.referenceExpression.toFixed(2)} CPM` : 'Unavailable'}
            {record.referenceExpressionStudies != null && record.referenceExpressionSamples != null && (
              <> · {record.referenceExpressionSamples} samples, {record.referenceExpressionStudies} studies</>
            )}
          </dd>
        </div>
        <div>
          <dt>Adjusted p-value</dt>
          <dd>{formatP(record.adjustedPValue)}</dd>
        </div>
        <div>
          <dt>Maximum change</dt>
          <dd>{record.valueAtMaxAbsChange.toFixed(2)}</dd>
        </div>
        <div>
          <dt>Shape cluster</dt>
          <dd>{record.cluster}</dd>
        </div>
        <div>
          <dt>Mean cells per donor</dt>
          <dd>{record.meanCellsPerDonor?.toFixed(1) ?? 'Unavailable'}</dd>
        </div>
      </dl>
      <div
        className="pair-card-annotations"
        aria-label={`Gene annotations from ${metadataRelease ?? 'PANTHER'}`}
      >
        <div>
          <span>Protein class</span>
          <p title={metadata?.proteinClasses.map((proteinClass) => proteinClass.label).join('; ')}>
            {metadata?.proteinClasses.length
              ? metadata.proteinClasses.map((proteinClass) => proteinClass.label).join('; ')
              : 'Missing'}
          </p>
        </div>
        <div>
          <span>PANTHER family</span>
          <p title={metadata?.families.map((family) => family.label).join('; ')}>
            {metadata?.families.length
              ? metadata.families.map((family) => family.label).join('; ')
              : 'Missing'}
          </p>
        </div>
        <small>Annotation: {metadataRelease ?? 'PANTHER metadata unavailable'}</small>
      </div>
      {(actions.onOpenGene || actions.onOpenCellType || actions.onOpenTrajectory) && (
        <div className="pair-card-actions">
          {actions.onOpenGene && (
            <button className="text-button" onClick={() => actions.onOpenGene?.(record)}>
              View gene profile →
            </button>
          )}
          {actions.onOpenCellType && (
            <button className="text-button" onClick={() => actions.onOpenCellType?.(record)}>
              View celltype profile →
            </button>
          )}
          {actions.onOpenTrajectory && (
            <button className="text-button" onClick={() => actions.onOpenTrajectory?.(record)}>
              Open cross-type trajectory map →
            </button>
          )}
        </div>
      )}
    </article>
  )
}
