const cpmFormatter = new Intl.NumberFormat('en', { maximumSignificantDigits: 4 })

export function formatCpm(value: number): string {
  return `${cpmFormatter.format(value)} CPM`
}
