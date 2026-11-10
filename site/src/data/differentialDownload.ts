import { AssetError, loadReleaseManifest } from './release'

export async function loadDifferentialDownload(path: string, releaseId: string, signal: AbortSignal): Promise<Blob> {
  const manifest = await loadReleaseManifest()
  if (manifest.releaseId !== releaseId)
    throw new AssetError('mismatch', 'The download belongs to a different release')
  if (!manifest.assets.some((asset) => asset.path === path && asset.kind === 'de_download'))
    throw new AssetError('unavailable', 'The download is not included in this release')
  const response = await fetch(
    `${import.meta.env.BASE_URL}data/${path}?release=${encodeURIComponent(releaseId)}`,
    { signal },
  )
  if (!response.ok || !response.body)
    throw new AssetError('unavailable', 'The download is unavailable')
  // Download assets are raw gzip files, served without HTTP Content-Encoding.
  // Legacy releases can still supply an ordinary CSV.
  const body = path.endsWith('.gz')
    ? response.body.pipeThrough(new DecompressionStream('gzip'))
    : response.body
  return new Response(body, { headers: { 'Content-Type': 'text/csv;charset=utf-8' } }).blob()
}
