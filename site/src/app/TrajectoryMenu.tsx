import { useEffect, useRef } from 'react'
import { prefetchReleaseAsset } from '../data/release'
import type { Page } from './navigation'

const pages: { id: Page; label: string; description: string }[] = [
  { id: 'genes', label: 'Genes', description: 'Compare a gene across cell types' },
  { id: 'cellTypes', label: 'Cell types', description: 'Explore genes within a cell type' },
  { id: 'trajectories', label: 'Cross-type trajectory map', description: 'Explore age-bin profile shapes' },
  { id: 'geneTypeMatrix', label: 'Gene × cell type', description: 'Compare age-bin changes for a gene set' },
]

export function TrajectoryMenu({ page, navigate }: { page: Page; navigate: (page: Page) => void }) {
  const menu = useRef<HTMLDetailsElement>(null)
  const active = pages.some((item) => item.id === page)
  useEffect(() => {
    const closeOutside = (event: PointerEvent) => {
      if (menu.current && !menu.current.contains(event.target as Node)) menu.current.open = false
    }
    document.addEventListener('pointerdown', closeOutside)
    return () => document.removeEventListener('pointerdown', closeOutside)
  }, [])
  return (
    <details
      ref={menu}
      className={`trajectory-menu ${active ? 'active' : ''}`}
      onPointerEnter={() => void prefetchReleaseAsset('pipeline-trajectory-map.json').catch(() => {})}
      onFocusCapture={() => void prefetchReleaseAsset('pipeline-trajectory-map.json').catch(() => {})}
      onBlur={(event) => {
        if (!event.currentTarget.contains(event.relatedTarget)) event.currentTarget.open = false
      }}
      onKeyDown={(event) => {
        if (event.key === 'Escape' && menu.current?.open) {
          event.preventDefault()
          menu.current.open = false
          menu.current.querySelector('summary')?.focus()
        }
      }}
    >
      <summary>
        Age-bin trajectories <span aria-hidden="true">⌄</span>
      </summary>
      <div className="trajectory-menu-panel">
        <p>EXPLORE AGE-BIN PROFILES</p>
        {pages.map((item) => (
          <button
            key={item.id}
            type="button"
            className={page === item.id ? 'active' : ''}
            aria-current={page === item.id ? 'page' : undefined}
            onClick={() => {
              if (menu.current) menu.current.open = false
              navigate(item.id)
            }}
          >
            <strong>{item.label}</strong>
            <span>{item.description}</span>
          </button>
        ))}
      </div>
    </details>
  )
}
