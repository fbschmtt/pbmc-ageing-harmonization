import { crc32, deflateSync, inflateSync } from 'node:zlib'
import { createHash } from 'node:crypto'
import { mkdir, writeFile } from 'node:fs/promises'
import { join } from 'node:path'

const signature = Buffer.from('89504e470d0a1a0a', 'hex')

// Recompress the existing filtered scanlines. Pixels, filters, alpha, colour
// profiles and all ancillary metadata stay unchanged. No image library needed.
export function optimizePng(bytes) {
  if (!bytes.subarray(0, 8).equals(signature)) throw new Error('Invalid PNG signature')
  const chunks = []
  for (let offset = 8; offset < bytes.length;) {
    if (offset + 12 > bytes.length) throw new Error('Truncated PNG chunk')
    const length = bytes.readUInt32BE(offset)
    if (offset + length + 12 > bytes.length) throw new Error('Truncated PNG payload')
    const type = bytes.toString('ascii', offset + 4, offset + 8)
    chunks.push({ type, bytes: bytes.subarray(offset, offset + length + 12), data: bytes.subarray(offset + 8, offset + 8 + length) })
    offset += length + 12
  }
  // Animated PNGs have separate frame streams and sequence numbers.
  if (chunks.some((chunk) => chunk.type === 'acTL')) return bytes
  const idat = chunks.filter((chunk) => chunk.type === 'IDAT')
  if (!idat.length || chunks.at(-1)?.type !== 'IEND') throw new Error('Incomplete PNG')
  const compressed = deflateSync(inflateSync(Buffer.concat(idat.map((chunk) => chunk.data))), { level: 9 })
  const replacement = Buffer.alloc(compressed.length + 12)
  replacement.writeUInt32BE(compressed.length)
  replacement.write('IDAT', 4)
  compressed.copy(replacement, 8)
  replacement.writeUInt32BE(crc32(replacement.subarray(4, -4)), replacement.length - 4)
  let inserted = false
  const result = Buffer.concat([signature, ...chunks.flatMap((chunk) => {
    if (chunk.type !== 'IDAT') return [chunk.bytes]
    if (inserted) return []
    inserted = true
    return [replacement]
  })])
  return result.length < bytes.length ? result : bytes
}

export async function writeSingleCellImage(stage, encoded) {
  const bytes = optimizePng(Buffer.from(encoded, 'base64'))
  const digest = createHash('sha256').update(bytes).digest('hex')
  const path = `single-cell-images/${digest}.png`
  await mkdir(join(stage, 'single-cell-images'), { recursive: true })
  await writeFile(join(stage, path), bytes)
  return path
}
