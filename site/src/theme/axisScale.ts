export type SymmetricAxis = {
  extent: number
  ticks: number[]
}

export type PositiveAxis = {
  maximum: number
  ticks: number[]
}

function tickStep(start: number, stop: number, targetIntervals: number) {
  const rawStep = Math.abs(stop - start) / Math.max(1, targetIntervals)
  const magnitude = 10 ** Math.floor(Math.log10(rawStep))
  const error = rawStep / magnitude
  const factor = error >= Math.sqrt(50) ? 10 : error >= Math.sqrt(10) ? 5 : error >= Math.sqrt(2) ? 2 : 1
  return factor * magnitude
}

function ticksBetween(start: number, stop: number, targetIntervals: number) {
  const step = tickStep(start, stop, targetIntervals)
  const first = Math.ceil(start / step)
  const last = Math.floor(stop / step)
  return Array.from({ length: Math.max(0, last - first + 1) }, (_, index) =>
    Number(((first + index) * step).toPrecision(12)),
  )
}

function finiteMaximum(values: number[]) {
  return Math.max(0, ...values.filter(Number.isFinite))
}

export function symmetricAxis(values: number[], paddingRatio = 0.05): SymmetricAxis {
  const largestMagnitude = finiteMaximum(values.map(Math.abs))
  const extent = largestMagnitude > 0 ? largestMagnitude * (1 + paddingRatio) : 1
  return { extent, ticks: ticksBetween(-extent, extent, 4) }
}

export function positiveAxis(values: number[], paddingRatio = 0.05): PositiveAxis {
  const largestValue = finiteMaximum(values)
  const maximum = largestValue > 0 ? largestValue * (1 + paddingRatio) : 1
  return { maximum, ticks: ticksBetween(0, maximum, 4) }
}

export function formatAxisValue(value: number, extent: number) {
  const normalized = Math.abs(value) < Math.abs(extent) * 1e-12 ? 0 : value
  const maximumFractionDigits = Math.min(8, Math.max(0, Math.ceil(-Math.log10(Math.abs(extent))) + 1))
  return new Intl.NumberFormat('en-US', {
    maximumFractionDigits,
    useGrouping: false,
  }).format(normalized)
}
