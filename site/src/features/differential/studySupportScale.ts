export function studySupportColor(significant: number, maximum: number): string {
  if (significant <= 0) return '#b9bdbc'
  const normalized = maximum > 0 ? Math.max(0, Math.min(1, significant / maximum)) : 0
  const lightness = 0.88 - normalized * 0.46
  return `oklch(${lightness.toFixed(3)} 0.09 180)`
}
