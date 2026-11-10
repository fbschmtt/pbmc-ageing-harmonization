import { useEffect, useRef, useState } from 'react'
import { loadDifferentialDownload } from '../../data/differentialDownload'

export function DifferentialDownload({ path, releaseId, filename }: {
  path: string
  releaseId: string
  filename: string
}) {
  const [busy, setBusy] = useState(false)
  const [failed, setFailed] = useState(false)
  const pending = useRef<AbortController | null>(null)
  useEffect(() => () => pending.current?.abort(), [])

  if (path.endsWith('.gz') && typeof DecompressionStream === 'undefined')
    return <>
      <a className="de-download-button"
        href={`${import.meta.env.BASE_URL}data/${path}?release=${encodeURIComponent(releaseId)}`}
        download={`${filename}.gz`}>Download compressed DEG results</a>
      <p>This browser cannot extract the download automatically. Extract the downloaded gzip archive to read the CSV.</p>
    </>

  const download = async () => {
    if (pending.current) return
    const controller = new AbortController()
    pending.current = controller
    setBusy(true)
    setFailed(false)
    try {
      const blob = await loadDifferentialDownload(path, releaseId, controller.signal)
      if (controller.signal.aborted) return
      const url = URL.createObjectURL(blob)
      const link = document.createElement('a')
      link.href = url
      link.download = filename
      document.body.append(link)
      link.click()
      link.remove()
      // Give the browser time to consume the URL before releasing its memory.
      window.setTimeout(() => URL.revokeObjectURL(url), 60_000)
    } catch {
      if (!controller.signal.aborted) setFailed(true)
    } finally {
      if (!controller.signal.aborted) {
        pending.current = null
        setBusy(false)
      }
    }
  }
  return <>
    <button className="de-download-button" type="button" disabled={busy} aria-busy={busy} onClick={download}>
      {busy ? 'Preparing CSV…' : 'Download DEG results'}
    </button>
    {busy && <p role="status">Preparing your CSV download…</p>}
    {failed && <p role="alert">Could not prepare the CSV download. Please try again.</p>}
  </>
}
