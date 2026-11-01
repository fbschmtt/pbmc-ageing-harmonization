import type { FractionAgeEstimate } from './FractionAgePlot'

export type StudyBaseline = { study: string; value: number; samples: number }

export function baselineSummary(estimates: FractionAgeEstimate[]) {
  const byStudy = new Map<string, StudyBaseline>()
  for (const item of estimates) {
    const value = item.referenceAdjustedFraction
    if (item.ageBin === item.referenceBin && value != null && Number.isFinite(value)) {
      byStudy.set(item.group, { study: item.group, value, samples: item.samples })
    }
  }
  const values = [...byStudy.values()].sort((a, b) => a.study.localeCompare(b.study))
  return values.length
    ? { mean: values.reduce((sum, item) => sum + item.value, 0) / values.length, studies: values }
    : null
}
