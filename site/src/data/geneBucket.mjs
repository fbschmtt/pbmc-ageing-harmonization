// Stable, case-insensitive partition for on-demand gene lookups.
export function geneBucket(gene) {
  let hash = 2166136261
  for (const character of gene.toUpperCase()) {
    hash ^= character.codePointAt(0)
    hash = Math.imul(hash, 16777619)
  }
  return (hash >>> 0) % 4
}
