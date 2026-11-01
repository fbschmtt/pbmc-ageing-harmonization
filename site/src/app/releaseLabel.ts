import type { ReleaseManifest } from '../data/release'

export function releaseLabel(release: ReleaseManifest | null): string {
  if (!release) return 'RELEASE UNAVAILABLE'
  const cells = release.inventory || []
  const synthetic =
    release.analysis?.syntheticTestData === true || cells.some((cell) => cell.syntheticTestData === true)
  const modes = new Set(cells.filter((cell) => cell.exported).map((cell) => cell.analysisMode || 'unknown'))
  if (release.analysis?.analysisMode) modes.add(release.analysis.analysisMode)
  const mode =
    modes.size > 1
      ? synthetic
        ? 'MIXED MODES · INCLUDES SYNTHETIC DATA'
        : 'MIXED ANALYSIS MODES'
      : synthetic
        ? 'SYNTHETIC TEST DATA'
        : modes.has('test_only')
          ? 'TEST DATA'
          : modes.has('standard')
            ? 'STANDARD ANALYSIS'
            : 'ANALYSIS MODE UNKNOWN'
  const date = /^\d{4}-\d{2}-\d{2}T/.test(release.generatedAt)
    ? release.generatedAt.slice(0, 10)
    : 'DATE UNKNOWN'
  return `${mode} · EXPORT ${date}`
}
