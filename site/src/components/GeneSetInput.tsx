import { useId, useRef, useState, type ReactNode } from 'react'
import type { PathwaySet } from '../data/pathways'
import './geneSetInput.css'

type GeneSetInputProps = {
  label: string
  value: string
  onChange: (value: string) => void
  helper: ReactNode
  className?: string
  pathwaySets?: PathwaySet[]
  activePathway?: PathwaySet | null
  onSelectPathway?: (pathway: PathwaySet) => void
  pathwayDisplay?: 'inline' | 'menu'
  onClear?: () => void
  clearLabel?: string
}

export function GeneSetInput({
  label,
  value,
  onChange,
  helper,
  className = '',
  pathwaySets = [],
  activePathway = null,
  onSelectPathway,
  pathwayDisplay = 'menu',
  onClear,
  clearLabel = 'Clear gene set',
}: GeneSetInputProps) {
  const [isPathwayMenuOpen, setIsPathwayMenuOpen] = useState(false)
  const menuId = useId()
  const triggerRef = useRef<HTMLButtonElement>(null)
  const closeMenu = () => {
    setIsPathwayMenuOpen(false)
    triggerRef.current?.focus()
  }
  const presets =
    pathwaySets.length > 0 && onSelectPathway ? (
      <div className={`gene-set-presets ${pathwayDisplay}`} aria-label="Predefined pathways">
        {pathwayDisplay === 'menu' ? (
          <div
            className="gene-set-menu"
            onBlur={(event) => {
              if (!event.currentTarget.contains(event.relatedTarget)) setIsPathwayMenuOpen(false)
            }}
            onKeyDown={(event) => {
              if (event.key === 'Escape') closeMenu()
            }}
          >
            <button
              ref={triggerRef}
              className="gene-set-menu-trigger"
              type="button"
              aria-expanded={isPathwayMenuOpen}
              aria-controls={menuId}
              onClick={() => setIsPathwayMenuOpen((open) => !open)}
            >
              <span>{activePathway?.label ?? 'Choose a predefined pathway'}</span>
              <span className="gene-set-menu-chevron" aria-hidden="true">
                ›
              </span>
            </button>
            <div
              id={menuId}
              className={`gene-set-preset-options ${isPathwayMenuOpen ? 'open' : ''}`}
              aria-hidden={!isPathwayMenuOpen}
            >
              {pathwaySets.map((pathway) => (
                <button
                  type="button"
                  key={pathway.id}
                  title={pathway.id}
                  tabIndex={isPathwayMenuOpen ? 0 : -1}
                  aria-pressed={activePathway?.id === pathway.id}
                  onClick={() => {
                    onSelectPathway(pathway)
                    closeMenu()
                  }}
                >
                  {pathway.label}
                </button>
              ))}
            </div>
          </div>
        ) : (
          pathwaySets.map((pathway) => (
            <button
              type="button"
              key={pathway.id}
              title={pathway.id}
              aria-pressed={activePathway?.id === pathway.id}
              onClick={() => onSelectPathway(pathway)}
              className={activePathway?.id === pathway.id ? 'active' : ''}
            >
              {pathway.label}
            </button>
          ))
        )}
      </div>
    ) : null
  return (
    <div className={`gene-set-input ${className}`}>
      <label>
        {label}
        {activePathway && <strong>{activePathway.label}</strong>}
        <textarea
          value={value}
          onChange={(event) => onChange(event.target.value)}
          placeholder="e.g. LAMP1, CTSD&#10;GBA"
          spellCheck={false}
        />
        <span className="gene-set-input-helper">{helper}</span>
      </label>
      {onClear && (
        <button className="gene-set-clear" type="button" disabled={!value.trim()} onClick={onClear}>
          {clearLabel}
        </button>
      )}
      {presets}
    </div>
  )
}
