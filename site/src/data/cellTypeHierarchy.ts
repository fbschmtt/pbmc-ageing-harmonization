import hierarchySource from '../../data-inputs/cell-type-hierarchy.json'

type CellTypeGroup = { label: string; cellTypes: string[] }

type MatchRule = {
  exact?: string[]
  prefixes?: string[]
  includes?: string[]
}

type HierarchyNode = {
  id: string
  label: string
  groupLabel?: string
  scope?: boolean
  match?: MatchRule
  children?: HierarchyNode[]
}

export type CompositionScope = {
  id: string
  label: string
  depth: number
  cellTypes: string[]
}

const hierarchy = hierarchySource as Record<string, HierarchyNode>

function matches(cellType: string, rule?: MatchRule) {
  return Boolean(
    rule?.exact?.includes(cellType) ||
      rule?.prefixes?.some((prefix) => cellType.startsWith(prefix)) ||
      rule?.includes?.some((fragment) => cellType.includes(fragment)),
  )
}

function leafNodes(node: HierarchyNode): HierarchyNode[] {
  if (!node.children?.length) return [node]
  return node.children.flatMap(leafNodes)
}

const leaves = Object.values(hierarchy).flatMap((tissue) => leafNodes(tissue))

function matchingLeaf(cellType: string) {
  return leaves.find((node) => matches(cellType, node.match))
}

export function groupCellTypes(cellTypes: string[]): CellTypeGroup[] {
  const groups = new Map<string, string[]>()
  for (const cellType of cellTypes) {
    const leaf = matchingLeaf(cellType)
    const label = leaf?.groupLabel ?? 'Other'
    groups.set(label, [...(groups.get(label) ?? []), cellType])
  }

  const orderedLabels = [
    ...leaves.map((node) => node.groupLabel).filter((label): label is string => Boolean(label)),
    'Other',
  ]
  return orderedLabels.flatMap((label) => {
    const items = groups.get(label)
    return items ? [{ label, cellTypes: items.sort() }] : []
  })
}

function cellTypesForNode(node: HierarchyNode, available: string[]): string[] {
  if (!node.children?.length) return available.filter((cellType) => matches(cellType, node.match)).sort()
  return node.children.flatMap((child) => cellTypesForNode(child, available))
}

function scopesForNode(node: HierarchyNode, available: string[], depth: number): CompositionScope[] {
  const assigned = cellTypesForNode(node, available)
  const cellTypes =
    depth === 0 ? [...assigned, ...available.filter((item) => !assigned.includes(item))] : assigned
  const current =
    node.scope === false || cellTypes.length === 0
      ? []
      : [{ id: node.id, label: node.label, depth, cellTypes }]
  return [
    ...current,
    ...(node.children?.flatMap((child) => scopesForNode(child, available, depth + 1)) ?? []),
  ]
}

export function compositionScopes(identities: { cellType: string; tissue: string }[]): CompositionScope[] {
  const byTissue = new Map<string, Set<string>>()
  for (const { cellType, tissue } of identities) {
    const values = byTissue.get(tissue) ?? new Set<string>()
    values.add(cellType)
    byTissue.set(tissue, values)
  }

  return Object.entries(hierarchy).flatMap(([tissue, node]) =>
    scopesForNode(node, [...(byTissue.get(tissue) ?? [])].sort(), 0),
  )
}
