import { useEffect, useRef } from 'react'
import type { NavPage, Page } from './navigation'

export function NavigationMenu({
  label,
  pages,
  page,
  navigate,
}: {
  label: string
  pages: NavPage[]
  page: Page
  navigate: (page: Page) => void
}) {
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
      className={`nav-menu ${active ? 'active' : ''}`}
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
        {label} <span aria-hidden="true">⌄</span>
      </summary>
      <div className="nav-menu-panel">
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
            {item.description && <span>{item.description}</span>}
          </button>
        ))}
      </div>
    </details>
  )
}
