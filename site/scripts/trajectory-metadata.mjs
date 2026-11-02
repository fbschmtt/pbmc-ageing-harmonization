// Older Python analysis outputs may contain non-standard non-finite literals.
// Preserve strings and let JSON.parse reject every other malformed input.
export function parseTrajectoryMetadata(text) {
  const normalized = text.replace(
    /"(?:\\.|[^"\\])*"|([:\[,]\s*)(?:NaN|Infinity|-Infinity)(?=\s*[,}\]])/g,
    (match, prefix) => (prefix === undefined ? match : `${prefix}null`),
  )
  return JSON.parse(normalized)
}
