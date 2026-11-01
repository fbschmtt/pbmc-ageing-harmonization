import type { FractionIndex } from '../../data/fractions'
import { baselineSummary } from './baselineSummary'
import { buildStudyColors } from '../../theme/studyColors'

export function fractionView(data: FractionIndex, cellType: string, denominator: string) {
  const records = data.records.filter((record) => record.cellType === cellType)
  const ageEstimates = data.ageEstimates.filter(
    (item) => item.cellType === cellType && item.denominator === denominator,
  )
  return {
    types: [...new Set(data.records.map((record) => record.cellType))].sort(),
    records,
    warning: records.find((record) => record.warning)?.warning,
    coefficients: data.coefficients.filter((item) => item.cellType === cellType),
    ageEstimates,
    skipped: data.ageStatus.filter(
      (item) =>
        item.cellType === cellType &&
        item.status === 'skipped' &&
        (!item.denominator || item.denominator === denominator),
    ),
    baseline: baselineSummary(ageEstimates),
    studyColors: buildStudyColors([
      ...data.records.map((item) => item.study),
      ...data.coefficients.map((item) => item.study),
      ...data.ageEstimates.map((item) => item.group),
    ]),
  }
}
