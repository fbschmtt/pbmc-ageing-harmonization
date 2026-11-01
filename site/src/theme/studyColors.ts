const studyPalette = [
  '#347a70',
  '#ca7559',
  '#6278a1',
  '#b39753',
  '#925b88',
  '#699154',
  '#a45559',
  '#447c9e',
  '#9d774e',
  '#704e95',
  '#607f70',
  '#9f6240',
  '#4b678c',
  '#8a7d4b',
]

/** Assign colors from the complete release study list, not the currently visible subset. */
export function buildStudyColors(studies: string[]): Map<string, string> {
  return new Map(
    [...new Set(studies)].sort().map((study, index) => [study, studyPalette[index % studyPalette.length]]),
  )
}
