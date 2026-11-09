const viridis = [
  '#440154',
  '#482878',
  '#3e4989',
  '#31688e',
  '#26828e',
  '#1f9e89',
  '#35b779',
  '#6ece58',
  '#b5de2b',
  '#fde725',
]

export function studySupportColor(significant: number, maximum: number): string {
  if (significant <= 0) return '#b9bdbc'
  const normalized = maximum > 1 ? Math.max(0, Math.min(1, (significant - 1) / (maximum - 1))) : 0
  return viridis[Math.round(normalized * (viridis.length - 1))]
}
