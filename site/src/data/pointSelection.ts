/** Prefer an available selection, otherwise choose the smallest finite adjusted p-value. */
export function choosePoint<T>(
  records: readonly T[],
  getPValue: (record: T) => number | null,
  preferred: (record: T) => boolean = () => false,
): T | null {
  function bestOf(pool: readonly T[]): T | null {
    let best: T | null = null
    let bestP = Infinity
    for (const record of pool) {
      const p = getPValue(record)
      if (p !== null && Number.isFinite(p) && p >= 0 && p <= 1 && p < bestP) {
        best = record
        bestP = p
      }
    }
    return best
  }
  return bestOf(records.filter(preferred)) ?? bestOf(records)
}
