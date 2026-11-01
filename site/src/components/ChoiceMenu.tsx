import { useEffect, useRef } from 'react'
import './choiceMenu.css'

export type ChoiceGroup = { label?: string; options: { value: string; label: string }[] }

export function ChoiceMenu({
  label,
  value,
  displayValue,
  groups,
  onChange,
}: {
  label: string
  value: string
  displayValue?: string
  groups: ChoiceGroup[]
  onChange: (value: string) => void
}) {
  const menu = useRef<HTMLDetailsElement>(null)
  const selected = groups.flatMap((group) => group.options).find((option) => option.value === value)
  useEffect(() => {
    const closeOutside = (event: PointerEvent) => {
      if (menu.current && !menu.current.contains(event.target as Node)) menu.current.open = false
    }
    document.addEventListener('pointerdown', closeOutside)
    return () => document.removeEventListener('pointerdown', closeOutside)
  }, [])
  return (
    <div className="choice-field">
      <span className="choice-label">{label}</span>
      <details
        className="choice-menu"
        ref={menu}
        onBlur={(event) => {
          if (!event.currentTarget.contains(event.relatedTarget)) event.currentTarget.open = false
        }}
        onKeyDown={(event) => {
          const current = menu.current
          if (!current) return
          if (event.key === 'Escape' && current.open) {
            event.preventDefault()
            current.open = false
            current.querySelector('summary')?.focus()
            return
          }
          if (!['ArrowDown', 'ArrowUp', 'Home', 'End'].includes(event.key)) return
          event.preventDefault()
          current.open = true
          const options = [...current.querySelectorAll<HTMLButtonElement>('[role="menuitemradio"]')]
          if (!options.length) return
          const index = options.indexOf(document.activeElement as HTMLButtonElement)
          const next =
            event.key === 'Home'
              ? 0
              : event.key === 'End'
                ? options.length - 1
                : event.key === 'ArrowDown'
                  ? (index + 1) % options.length
                  : index <= 0
                    ? options.length - 1
                    : index - 1
          options[next].focus()
        }}
      >
        <summary
          aria-label={`${label}: ${displayValue || selected?.label || 'Select an option'}`}
          aria-haspopup="menu"
        >
          <span className="choice-value">{displayValue || selected?.label || 'Select an option'}</span>
          <span className="choice-chevron" aria-hidden="true">
            ⌄
          </span>
        </summary>
        <div className="choice-panel" role="menu" aria-label={label}>
          {groups.map((group, groupIndex) => (
            <div key={`${group.label || ''}-${groupIndex}`} role="group" aria-label={group.label}>
              {group.label && <p className="choice-group">{group.label}</p>}
              {group.options.map((option, index) => (
                <button
                  key={`${option.value}-${index}`}
                  type="button"
                  role="menuitemradio"
                  aria-checked={option.value === value}
                  onClick={() => {
                    onChange(option.value)
                    if (menu.current) {
                      menu.current.open = false
                      menu.current.querySelector('summary')?.focus()
                    }
                  }}
                >
                  <span>{option.label}</span>
                  {option.value === value && (
                    <span className="choice-check" aria-hidden="true">
                      ✓
                    </span>
                  )}
                </button>
              ))}
            </div>
          ))}
        </div>
      </details>
    </div>
  )
}
