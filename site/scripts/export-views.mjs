import { exportTrajectoryView } from './export-trajectory-view.mjs'
import { exportFractionView } from './export-fraction-view.mjs'

export async function exportViews(input, output, cellTypes) {
  const counts = {
    trajectory: { input: 0, exported: 0 },
    fractions: { input: 0, exported: 0 },
    coefficients: { input: 0, exported: 0 },
    ageEstimates: { input: 0, exported: 0 },
  }
  console.log('[site-export] Preparing trajectory and reference-expression views')
  await exportTrajectoryView(input, output, cellTypes, counts)
  console.log('[site-export] Preparing cell-fraction and age-bin summaries')
  await exportFractionView(input, output, cellTypes, counts)
  console.log(
    `Exported ${counts.trajectory.exported} trajectory records and ${counts.fractions.exported} study fraction summaries`,
  )
  return counts
}
