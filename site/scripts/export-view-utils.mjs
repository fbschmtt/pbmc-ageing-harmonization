import { readFile, unlink } from 'node:fs/promises'
export const num = (value) =>
  value == null || value.trim() === '' || !Number.isFinite(Number(value)) ? null : Number(value)
export async function optional(path) {
  try {
    return await readFile(path, 'utf8')
  } catch (error) {
    if (error.code === 'ENOENT') return null
    throw error
  }
}
export async function remove(path) {
  try {
    await unlink(path)
  } catch (error) {
    if (error.code !== 'ENOENT') throw error
  }
}

export function exclude(count, reason, amount = 1) {
  count.excluded ??= {}
  count.excluded[reason] = (count.excluded[reason] || 0) + amount
}
