import { useEffect, useState } from 'react'
import { navigation, type Page } from './navigation'

function pageFromHash(): Page {
  const page = window.location.hash.slice(1) as Page
  return navigation.some((item) => item.id === page) ? page : 'overview'
}

export function usePageNavigation() {
  const [page, setPage] = useState<Page>(pageFromHash)
  useEffect(() => {
    const syncPage = () => setPage(pageFromHash())
    window.addEventListener('hashchange', syncPage)
    if (!window.location.hash) window.history.replaceState(null, '', '#overview')
    return () => window.removeEventListener('hashchange', syncPage)
  }, [])
  const navigate = (nextPage: Page) => {
    window.location.hash = nextPage
    setPage(nextPage)
  }
  return { page, navigate }
}
