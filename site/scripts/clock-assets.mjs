function assetSegment(value) {
  const segment = value.replace(/[^A-Za-z0-9._-]/g, '_')
  if (!segment || segment === '.' || segment === '..') throw new Error(`Invalid clock asset identifier: ${value}`)
  return segment
}

export function clockModelAssetPath(cellTypeId, methodId) {
  return `clock-models/${assetSegment(cellTypeId)}/${assetSegment(methodId)}.json`
}
