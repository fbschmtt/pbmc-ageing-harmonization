// AIFI Immune Health Atlas L2 palette: https://apps.allenimmunology.org/aifi/resources/imm-health-atlas/downloads/models/
const colors: Record<string, string> = {
  'Naive CD4 T cell': '#597FC6',
  'Memory CD4 T cell': '#233C5B',
  'Treg': '#001A9A',
  'DN T cell': '#656E68',
  'Proliferating T cell': '#AEAEB0',
  'Naive CD8 T cell': '#63A686',
  'Memory CD8 T cell': '#314937',
  'CD8aa': '#7ECDAD',
  'gdT': '#72AF5A',
  'MAIT': '#1E8F64',
  'Transitional B cell': '#9C6469',
  'Naive B cell': '#D5ABAB',
  'Memory B cell': '#B16B72',
  'Effector B cell': '#F4CEC3',
  'Plasma cell': '#784D47',
  'CD56bright NK cell': '#FF7F55',
  'CD56dim NK cell': '#E16040',
  'Proliferating NK cell': '#CC3B18',
  'CD14 monocyte': '#F6831B',
  'Intermediate monocyte': '#B75228',
  'CD16 monocyte': '#FCAF87',
  'cDC2': '#BE6E23',
  'cDC1': '#F8C755',
  'pDC': '#FBB64E',
  'ASDC': '#8A4E1C',
  'ILC': '#F4AD3A',
  'Progenitor cell': '#A31001',
  'Platelet': '#5B3930',
  'Erythrocyte': '#1F1515',
}

export function cellTypeColor(label: string): string | undefined {
  return colors[label]
}
