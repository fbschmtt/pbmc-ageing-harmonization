const viridisStops = [
  [68, 1, 84],
  [71, 44, 122],
  [59, 81, 139],
  [44, 113, 142],
  [33, 144, 141],
  [39, 173, 129],
  [92, 200, 99],
  [170, 220, 50],
  [253, 231, 37],
] as const

export function studySupportColor(significant: number, maximum: number): string {
  const normalized = maximum > 0 ? Math.max(0, Math.min(1, significant / maximum)) : 0
  const scaled = normalized * (viridisStops.length - 1)
  const lower = Math.floor(scaled)
  const upper = Math.min(viridisStops.length - 1, lower + 1)
  const fraction = scaled - lower
  const channels = viridisStops[lower].map((value, index) =>
    Math.round(value + (viridisStops[upper][index] - value) * fraction),
  )
  return `rgb(${channels.join(' ')})`
}
