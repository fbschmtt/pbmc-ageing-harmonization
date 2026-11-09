import { useEffect, useState } from 'react'
import { isPage, type Page } from './navigation'

function pageFromHash(): Page {
  const page = window.location.hash.slice(1) as Page
  return isPage(page) ? page : 'overview'
}

export function usePageNavigation() {
  const [page, setPage] = useState<Page>(pageFromHash)
  useEffect(() => {
    const syncPage = () => setPage(pageFromHash())
    window.addEventListener('hashchange', syncPage)
    if (!window.location.hash) window.history.replaceState(null, '', '#overview')
    return () => window.removeEventListener('hashchange', syncPage)
  }, [])
  useEffect(() => {
    window.scrollTo(0, 0)
  }, [page])
  const navigate = (nextPage: Page) => {
    window.location.hash = nextPage
    setPage(nextPage)
  }
  return { page, navigate }
}
