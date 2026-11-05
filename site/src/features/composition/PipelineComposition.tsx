import { cellTypeColor } from '../../theme/cellTypeColors'
import { assetFailureKind, type AssetFailureKind } from '../../data/release'
import { loadFractions, type FractionIndex } from '../../data/fractions'
import { useEffect, useState } from 'react'
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
        {failure && <p>Study-specific fraction summaries could not be loaded.</p>}
      </section>
    )
  const hasPreferred =
    preferredCellType !== null && data.records.some((record) => record.cellType === preferredCellType)
  const selectedCellType = hasPreferred ? preferredCellType : cellType
  const { types, records, warning, coefficients, ageEstimates, skipped, studyColors } =
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
    </section>
  )
}
