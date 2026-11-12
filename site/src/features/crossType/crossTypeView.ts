import type { CrossTypeIndex, CrossTypeGene } from '../../data/crossType'

export function visibleGenes(index: CrossTypeIndex, query: string, minimum: number) {
  const text = query.trim().toUpperCase()
  return index.genes.filter((row) => row.significant >= minimum && row.gene.toUpperCase().includes(text))
}
export function orderedResults(index: CrossTypeIndex, detail: CrossTypeGene) {
  const records = new Map(detail.records.map((row) => [row.cellType, row]))
  return index.cellTypes.map((cellType) => ({ cellType, result: records.get(cellType) ?? null }))
}
export function patternExtent(patterns: { mean: number[] }[]) {
  return Math.max(0.25, ...patterns.flatMap((pattern) => pattern.mean.map(Math.abs)))
}
export function patternPoints(values: number[], maximum = Math.max(0.25, ...values.map(Math.abs))) {
  return values
    .map((value, i) => `${30 + (i / Math.max(1, values.length - 1)) * 260},${78 - (value / maximum) * 50}`)
    .join(' ')
}
