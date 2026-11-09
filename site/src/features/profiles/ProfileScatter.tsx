import { cellTypeColor } from '../../theme/cellTypeColors'
import { VolcanoPlot } from '../../components/VolcanoPlot'
import type { TrajectoryRecord } from '../../data/atlas'
import { visualColors } from '../../theme/visualEncoding'
import { formatNumber } from '../../theme/format'

type ProfileScatterProps = {
  records: TrajectoryRecord[]
  selected: TrajectoryRecord | null
  onSelect: (record: TrajectoryRecord) => void
  label: (record: TrajectoryRecord) => string
  colorByCellType?: boolean
  highlightedGenes?: ReadonlySet<string>
}
export function ProfileScatter({
  records,
  selected,
  onSelect,
  label,
  highlightedGenes,
  colorByCellType = false,
}: ProfileScatterProps) {
  const highlighted = (record: TrajectoryRecord) => highlightedGenes?.has(record.gene.toUpperCase()) ?? false
  return (
    <VolcanoPlot
      className="profile-scatter"
      records={records}
      selected={selected}
      onSelect={onSelect}
      id={(record) => record.id}
      effect={(record) => record.valueAtMaxAbsChange}
      adjustedP={(record) => record.adjustedPValue}
      highlighted={highlighted}
      xLabel="Signed maximum age-bin change (log₂ fold change)"
      ariaLabel="Signed maximum age-bin change against negative log10 adjusted p-value"
      title={(record) =>
        `${label(record)} · ${record.tissue}\nSigned maximum age change: ${formatNumber(record.valueAtMaxAbsChange)}\nAdjusted p-value: ${formatNumber(record.adjustedPValue)}`
      }
      appearance={(record, isSelected) => ({
        fill: colorByCellType ? (cellTypeColor(record.cellType) ?? '#87918d') : record.tissue === 'brain' ? visualColors.brain : visualColors.blood,
        stroke: isSelected
          ? visualColors.selected
          : highlighted(record)
            ? visualColors.geneSetOutline
            : 'none',
        strokeWidth: isSelected ? 2.5 : highlighted(record) ? 1.5 : 0,
        radius: isSelected ? 7.5 : highlighted(record) ? 6 : 5,
      })}
    />
  )
}
