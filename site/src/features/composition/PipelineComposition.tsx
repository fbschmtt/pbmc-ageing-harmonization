import { cellTypeColor } from '../../theme/cellTypeColors'
import { assetFailureKind, type AssetFailureKind } from '../../data/release'
import { loadFractions, type FractionIndex } from '../../data/fractions'
import { useEffect, useState } from 'react'
import { BaselinePlot } from './BaselinePlot'
import { fractionView } from './fractionView'
import { FractionAgePlot } from './FractionAgePlot'
import { FractionForest } from './FractionForest'
import { ChoiceMenu } from '../../components/ChoiceMenu'
import { groupCellTypes } from '../../data/cellTypeHierarchy'
import './composition.css'

export function PipelineComposition({
  preferredCellType,
  onSelectCellType,
  denominator,
  onSelectDenominator,
}: {
  preferredCellType: string | null
  onSelectCellType: (cellType: string) => void
  denominator: string
  onSelectDenominator: (denominator: string) => void
}) {
  const [data, setData] = useState<FractionIndex | null>(null)
  const [failure, setFailure] = useState<AssetFailureKind | null>(null)
  const [cellType, setCellType] = useState('')
  useEffect(() => {
    let active = true
    loadFractions()
      .then((payload) => {
        if (!active) return
        setData(payload)
        setCellType(payload.records[0].cellType)
      })
      .catch((error) => {
        if (active) setFailure(assetFailureKind(error))
      })
    return () => {
      active = false
    }
  }, [])
  if (!data)
    return (
      <section className="content-page">
        <p className="eyebrow">PBMC composition</p>
        <h1>
          {failure === 'invalid'
            ? 'Fraction data is invalid.'
            : failure === 'mismatch'
              ? 'Fraction data belongs to a different release.'
              : failure
                ? 'Fraction summaries are not available.'
                : 'Loading fraction summaries…'}
        </h1>
        <p>
          Export results from completed per-type DE reports with <code>make site-export</code>. The reports
          provide the study-specific fraction diagnostics used here.
        </p>
      </section>
    )
  const hasPreferred =
    preferredCellType !== null && data.records.some((record) => record.cellType === preferredCellType)
  const selectedCellType = hasPreferred ? preferredCellType : cellType
  const { types, records, warning, coefficients, ageEstimates, skipped, baseline, studyColors } =
    fractionView(data, selectedCellType, denominator)
  return (
    <section className="content-page pipeline-composition">
      <h1 className="visually-hidden">Cell fractions, by study</h1>
      {warning && (
        <div className="de-warning" role="status">
          {warning}
        </div>
      )}
      {preferredCellType && !hasPreferred && (
        <p role="status">
          Fraction summaries are unavailable for {preferredCellType}. Showing {selectedCellType}.{' '}
          <button type="button" onClick={() => onSelectCellType(selectedCellType)}>
            Use {selectedCellType}
          </button>
        </p>
      )}
      <ChoiceMenu
        label="Cell type"
        value={selectedCellType}
        groups={groupCellTypes(types).map((group) => ({
          label: group.label,
          options: group.cellTypes.map((type) => ({ value: type, label: type, color: cellTypeColor(type) })),
        }))}
        onChange={(value) => {
          setCellType(value)
          onSelectCellType(value)
        }}
      />
      <FractionForest coefficients={coefficients} diagnostics={records} studyColors={studyColors} />
      <ChoiceMenu
        label="Age plot denominator"
        value={denominator}
        groups={[
          {
            options: [
              { value: 'l1_parent', label: 'Within parent compartment' },
              { value: 'retained_pbmc', label: 'Within retained PBMCs' },
            ],
          },
        ]}
        onChange={onSelectDenominator}
      />
      <FractionAgePlot estimates={ageEstimates} studyColors={studyColors} />
      {skipped.map((item, index) => (
        <p key={`${item.study}-${index}`}>
          {item.study || 'Age-bin model'}: {item.reason}
        </p>
      ))}
      <div className="de-table-wrap">
        <table>
          <thead>
            <tr>
              <th>Study</th>
              <th>Parent denominator</th>
              <th>Mean fraction</th>
              <th>Adult samples</th>
              <th>Model covariates</th>
            </tr>
          </thead>
          <tbody>
            {records.map((record) => (
              <tr key={`${record.cellType}-${record.study}`}>
                <td>{record.study}</td>
                <td>{record.parent}</td>
                <td>{(record.meanFraction * 100).toFixed(2)}%</td>
                <td>{record.samples}</td>
                <td>{record.covariates?.replaceAll(';', ', ') || '—'}</td>
              </tr>
            ))}
          </tbody>
        </table>
      </div>
      <p>
        Fractions describe relative representation within the captured parent compartment; they do not measure
        absolute cell abundance. Studies are shown separately, and may include different covariates and
        samples. Effects are shown in percentage points, with age measured per decade.
      </p>
      <div className="fraction-baseline" role="status">
        <span className="eyebrow">Model-adjusted reference fraction</span>
        {baseline ? (
          <>
            <strong>{(baseline.mean * 100).toFixed(2)}%</strong>
            <p>
              Equal-weight mean across {baseline.studies.length}{' '}
              {baseline.studies.length === 1 ? 'study' : 'studies'} in the {ageEstimates[0].referenceBin}{' '}
              reference group. Each study's value averages fitted predictions across its model samples,
              preserving their mean BMI and the proportions of included categorical covariates.
            </p>
            <BaselinePlot studies={baseline.studies} mean={baseline.mean} studyColors={studyColors} />
            <p>
              The dashed line is the equal-weight study mean. Dots are study estimates, not individual
              samples; uncertainty intervals are not shown. Linear-model predictions can fall outside 0–100%
              and are not clipped.
            </p>
          </>
        ) : (
          <p>
            Adjusted reference fractions are unavailable. Rerender the DE reports and export the site data.
          </p>
        )}
      </div>
    </section>
  )
}
