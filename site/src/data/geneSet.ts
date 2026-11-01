export function parseGeneSet(value: string): Set<string> {
  const genes = new Set<string>()
  let gene = ''
  const addGene = () => {
    if (gene) genes.add(gene.toUpperCase())
    gene = ''
  }
  for (const character of value) {
    if (character === ',' || character.trim() === '') addGene()
    else gene += character
  }
  addGene()
  return genes
}
