import { expect, it } from 'vitest'
import { releaseLabel } from './releaseLabel'
import type { ReleaseManifest } from '../data/release'

const base: ReleaseManifest = {
  schemaVersion: 1,
  releaseId: 'export-1',
  generatedAt: '2026-10-31T12:00:00Z',
  assets: [],
}
it('labels synthetic, standard, mixed and unknown modes without inventing an analysis version', () => {
  expect(releaseLabel({ ...base, analysis: { syntheticTestData: true, analysisMode: 'test_only' } })).toBe(
    'SYNTHETIC TEST DATA · EXPORT 2026-10-31',
  )
  expect(releaseLabel({ ...base, analysis: { analysisMode: 'standard' } })).toBe(
    'STANDARD ANALYSIS · EXPORT 2026-10-31',
  )
  expect(
    releaseLabel({
      ...base,
      inventory: [
        { id: 'a', label: 'A', analysisMode: 'standard', exported: true },
        { id: 'b', label: 'B', analysisMode: 'test_only', exported: true },
      ],
    }),
  ).toBe('MIXED ANALYSIS MODES · EXPORT 2026-10-31')
  expect(releaseLabel(base)).toBe('ANALYSIS MODE UNKNOWN · EXPORT 2026-10-31')
})
