import { cellTypeColor } from '../theme/cellTypeColors'
import { ChoiceMenu } from './ChoiceMenu'

export type CellTypeOption = { value: string; label: string; color?: string }
export type CellTypeOptionGroup = { label?: string; options: CellTypeOption[] }

export function CellTypeSelector({
  label = 'Cell type',
  value,
  groups,
  onChange,
  className,
}: {
  label?: string
  value: string
  groups: CellTypeOptionGroup[]
  onChange: (value: string) => void
  className?: string
}) {
  return (
    <div className={className}>
      <ChoiceMenu
        label={label}
        value={value}
        groups={groups.map((group) => ({
          ...group,
          options: group.options.map((option) => ({
            ...option,
            color: option.color ?? cellTypeColor(option.label),
          })),
        }))}
        onChange={onChange}
      />
    </div>
  )
}
