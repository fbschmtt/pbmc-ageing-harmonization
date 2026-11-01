export const goTermLimits = {
  biologicalProcesses: 3,
  molecularFunctions: 2,
  cellularComponents: 2,
}

const uninformativeGoIds = new Set([
  'GO:0003674', // molecular function
  'GO:0008150', // biological process
  'GO:0005575', // cellular component
  'GO:0110165', // cellular anatomical entity
  'GO:0005622', // intracellular anatomical structure
])

export function selectInformativeTerms(terms, frequencies, limit) {
  return [...terms.values()]
    .filter((term) => !uninformativeGoIds.has(term.id))
    .sort(
      (left, right) =>
        (frequencies.get(left.id) ?? Number.POSITIVE_INFINITY) -
          (frequencies.get(right.id) ?? Number.POSITIVE_INFINITY) || left.label.localeCompare(right.label),
    )
    .slice(0, limit)
}
