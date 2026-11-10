import { crc32, deflateSync, inflateSync } from 'node:zlib'
import { expect, it } from 'vitest'
import { optimizePng } from './png.mjs'

function chunk(type, data) {
  const result = Buffer.alloc(data.length + 12)
  result.writeUInt32BE(data.length)
  result.write(type, 4)
  data.copy(result, 8)
  result.writeUInt32BE(crc32(result.subarray(4, -4)), result.length - 4)
  return result
}
function fixture(level = 0, split = true) {
  const signature = Buffer.from('89504e470d0a1a0a', 'hex')
  const header = Buffer.alloc(13)
  header.writeUInt32BE(128, 0); header.writeUInt32BE(128, 4)
  header[8] = 8; header[9] = 6 // RGBA; retain alpha and bit depth.
  const pixels = Buffer.alloc(128 * (1 + 128 * 4), 0)
  const compressed = deflateSync(pixels, { level })
  const midpoint = split ? Math.floor(compressed.length / 2) : compressed.length
  const metadata = [chunk('tEXt', Buffer.from('Description\0Scientific figure')), chunk('pHYs', Buffer.from('00000ec400000ec401', 'hex'))]
  const bytes = Buffer.concat([signature, chunk('IHDR', header), ...metadata,
    chunk('IDAT', compressed.subarray(0, midpoint)),
    ...(split ? [chunk('IDAT', compressed.subarray(midpoint))] : []), chunk('IEND', Buffer.alloc(0))])
  return { bytes, pixels }
}
function readChunks(bytes) {
  const chunks = []
  for (let offset = 8; offset < bytes.length;) {
    const length = bytes.readUInt32BE(offset)
    const raw = bytes.subarray(offset, offset + 12 + length)
    expect(raw.readUInt32BE(raw.length - 4)).toBe(crc32(raw.subarray(4, -4)))
    chunks.push({ type: raw.toString('ascii', 4, 8), data: raw.subarray(8, -4), raw })
    offset += length + 12
  }
  return chunks
}
it('automatically reduces PNG storage while preserving exact scanlines and all other chunks', () => {
  const { bytes, pixels } = fixture()
  const optimized = optimizePng(bytes)
  expect(optimized.length).toBeLessThan(bytes.length)
  const before = readChunks(bytes), after = readChunks(optimized)
  expect(inflateSync(Buffer.concat(after.filter((c) => c.type === 'IDAT').map((c) => c.data)))).toEqual(pixels)
  expect(after.filter((c) => c.type !== 'IDAT')).toEqual(before.filter((c) => c.type !== 'IDAT'))
})
it('keeps the original file when the stronger pass gives no saving', () => {
  const { bytes } = fixture(9, false)
  expect(optimizePng(bytes)).toBe(bytes)
})
it('rejects truncated PNG data instead of publishing a partial figure', () => {
  expect(() => optimizePng(fixture().bytes.subarray(0, -5))).toThrow('Truncated')
})
