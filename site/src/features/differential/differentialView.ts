import type { Cell, Gene, Model } from '../../data/differential'
import { defaultModelId } from './modelSelection'

export function hasPlottableEstimate(record: Gene): record is Gene & { effect: number; padj: number } {
  return (
    record.effect !== null &&
    Number.isFinite(record.effect) &&
    record.padj !== null &&
    Number.isFinite(record.padj) &&
    record.padj >= 0 &&
    record.padj <= 1
  )
}

export function visibleGenes(records: Gene[], requested: ReadonlySet<string>): Gene[] {
  return records.filter(
    (record) => (record.padj !== null && record.padj < 0.05) || requested.has(record.gene.toUpperCase()),
  )
}

export function initialCell(cells: Cell[], gene: string | null, type: string | null): Cell {
  const hasGene = (item: Cell) =>
    item.models
      .find((model) => model.id === defaultModelId(item.models))
      ?.records.some((record) => record.gene === gene && hasPlottableEstimate(record))
  return cells.find((item) => item.label === type) ?? cells.find(hasGene) ?? cells[0]
}

export function filterGenes(model: Model | undefined, search: string): Gene[] {
  const query = search.trim().toUpperCase()
  return (model?.records || [])
    .filter(hasPlottableEstimate)
    .filter((record) => record.gene.toUpperCase().includes(query))
    .sort((a, b) => (a.padj ?? Infinity) - (b.padj ?? Infinity) || a.gene.localeCompare(b.gene))
}
