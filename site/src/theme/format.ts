export function formatNumber(value: number | null | undefined): string {
  if (value == null || !Number.isFinite(value)) return '—'
  if (value === 0) return '0'
  return Math.abs(value) >= 1 ? value.toFixed(2) : value.toPrecision(2)
}

export function formatCpm(value: number): string {
  return `${formatNumber(value)} CPM`
}
