import { cellTypeColor } from '../../theme/cellTypeColors'
import { VolcanoPlot } from '../../components/VolcanoPlot'
import type { Gene } from '../../data/contracts.mjs'
import { visualColors } from '../../theme/visualEncoding'
import { formatNumber } from '../../theme/format'

export type CombinedAgePoint = Gene & {
  id: string
  cellType: string
  cellId: string
  effect: number
  padj: number
  trajectoryAsset: string | null
}
export function ProfileScatter({ records, selected, onSelect }: {
  records: CombinedAgePoint[]
  selected: CombinedAgePoint | null
  onSelect: (record: CombinedAgePoint) => void
}) {
  return <VolcanoPlot className="profile-scatter" aspectRatio="square" threshold={0.05} records={records} selected={selected} onSelect={onSelect}
    id={(row) => row.id} effect={(row) => row.effect * 10} adjustedP={(row) => row.padj}
    xLabel="Combined linear age effect (log₂ fold change per decade)"
    ariaLabel="Combined linear age effect per decade against negative log10 adjusted p-value"
    title={(row) => `${row.cellType} · pbmc\nCombined linear age effect per decade: ${formatNumber(row.effect * 10)}\nAdjusted p-value (combined age model): ${formatNumber(row.padj)}`}
    appearance={(row, isSelected) => ({ fill: cellTypeColor(row.cellType) ?? '#87918d',
      stroke: isSelected ? visualColors.selected : 'none', strokeWidth: isSelected ? 2.5 : 0, radius: isSelected ? 7.5 : 5 })} />
}
