import { describe, expect, it } from 'vitest'
import type { TrajectoryRecord } from '../data/atlas'
import type { PathwaySet } from '../data/pathways'
import { atlasSessionReducer, initialAtlasSession } from './session'

const pathway = (id: string, genes: string[]): PathwaySet => ({
  id,
  genes,
  label: id,
  description: id,
  sourceFile: `${id}.gmt`,
})

const record: TrajectoryRecord = {
  id: 'HK1-Astro',
  gene: 'HK1',
  cellType: 'Astro',
  tissue: 'brain',
  x: 1,
  y: 2,
  adjustedPValue: 0.01,
  referenceExpression: 1,
  meanCellsPerDonor: 10,
  valueAtMaxAbsChange: 0.5,
  cluster: '1',
  trajectory: [0, 0, 0, 0, 0, 0, 0],
}

describe('atlasSessionReducer', () => {
  it('retains the selected gene when switching other session controls', () => {
    const selected = atlasSessionReducer(initialAtlasSession, { type: 'selectGene', gene: 'HK1' })
    const edited = atlasSessionReducer(selected, { type: 'editGeneSet', value: 'GBA' })
    expect(edited.selectedGene).toBe('HK1')
  })
  it('updates gene and cell context and clears stale selections in other views', () => {
    const previous = {
      ...initialAtlasSession,
      selectedGene: 'OLD',
      selectedCellType: 'Old type',
      geneProfile: { gene: 'OLD', cellType: 'Old type', selectedId: 'old-gene' },
      cellTypeProfile: { gene: 'OLD', cellType: 'Old type', selectedId: 'old-type' },
      trajectory: { ...initialAtlasSession.trajectory, query: 'OLD', selectedId: 'old-map' },
    }
    const selected = atlasSessionReducer(previous, { type: 'selectGene', gene: 'HK1', cellType: 'Astro' })
    expect(selected.selectedGene).toBe('HK1')
    expect(selected.selectedCellType).toBe('Astro')
    expect(selected.geneProfile).toMatchObject({ gene: 'HK1', selectedId: null })
    expect(selected.cellTypeProfile).toMatchObject({ cellType: 'Astro', selectedId: null })
    expect(selected.trajectory).toMatchObject({ query: '', selectedId: null })
  })

  it('shares a cell-type selection without losing the selected gene', () => {
    const previous = {
      ...initialAtlasSession,
      selectedGene: 'HK1',
      selectedCellType: 'Old',
      geneProfile: { gene: 'HK1', cellType: 'Old', selectedId: 'old-pair' },
      cellTypeProfile: { gene: 'HK1', cellType: 'Old', selectedId: 'old-pair' },
      trajectory: { ...initialAtlasSession.trajectory, query: 'HK1', selectedId: 'old-pair' },
    }
    const selected = atlasSessionReducer(previous, { type: 'selectCellType', cellType: 'Astro' })
    expect(selected.selectedGene).toBe('HK1')
    expect(selected.selectedCellType).toBe('Astro')
    expect(selected.geneProfile.selectedId).toBeNull()
    expect(selected.cellTypeProfile).toMatchObject({ cellType: 'Astro', selectedId: null })
    expect(selected.trajectory).toMatchObject({ query: 'HK1', selectedId: null })
    expect(atlasSessionReducer(selected, { type: 'selectCellType', cellType: 'Astro' })).toBe(selected)
  })

  it('keeps the model and composition denominator when changing other controls', () => {
    const withModel = atlasSessionReducer(initialAtlasSession, {
      type: 'selectDEModel',
      modelId: 'per_study:study_a:age',
    })
    const withDenominator = atlasSessionReducer(withModel, {
      type: 'selectCompositionDenominator',
      denominator: 'retained_pbmc',
    })
    const selected = atlasSessionReducer(withDenominator, { type: 'selectCellType', cellType: 'Astro' })
    expect(selected).toMatchObject({
      selectedDEModelId: 'per_study:study_a:age',
      compositionDenominator: 'retained_pbmc',
      selectedCellType: 'Astro',
    })
  })

  it('prefers Cellular senescence for first-load initialization', () => {
    const state = atlasSessionReducer(initialAtlasSession, {
      type: 'initializePathway',
      pathways: [pathway('FIRST', ['A']), pathway('REACTOME_CELLULAR_SENESCENCE', ['B', 'C'])],
    })
    expect(state.activePathwayId).toBe('REACTOME_CELLULAR_SENESCENCE')
    expect(state.highlightGenes).toBe('B\nC')
  })

  it('does not overwrite a manual edit when pathway loading finishes later', () => {
    const edited = atlasSessionReducer(initialAtlasSession, { type: 'editGeneSet', value: 'LAMP1' })
    const initialized = atlasSessionReducer(edited, {
      type: 'initializePathway',
      pathways: [pathway('REACTOME_CELLULAR_SENESCENCE', ['B'])],
    })
    expect(initialized).toBe(edited)
  })

  it('tracks a selected preset until the shared gene set is edited', () => {
    const selected = atlasSessionReducer(initialAtlasSession, {
      type: 'selectPathway',
      pathway: pathway('MITOPHAGY', ['PINK1', 'PRKN']),
    })
    expect(selected).toMatchObject({ activePathwayId: 'MITOPHAGY', highlightGenes: 'PINK1\nPRKN' })
    expect(atlasSessionReducer(selected, { type: 'editGeneSet', value: 'PINK1' })).toMatchObject({
      activePathwayId: null,
      highlightGenes: 'PINK1',
    })
  })

  it('keeps gene and cell-type profile memory independent', () => {
    const geneState = atlasSessionReducer(initialAtlasSession, {
      type: 'openProfile',
      tab: 'gene',
      record,
      selectPair: true,
    })
    const cellTypeState = atlasSessionReducer(geneState, {
      type: 'openProfile',
      tab: 'cellType',
      record: { ...record, id: 'GBA-Micro', gene: 'GBA', cellType: 'Micro' },
      selectPair: false,
    })
    expect(cellTypeState.geneProfile.selectedId).toBe('HK1-Astro')
    expect(cellTypeState.cellTypeProfile).toMatchObject({ cellType: 'Micro', selectedId: null })
  })

  it('opens a requested trajectory with a clean search while preserving the UMAP view', () => {
    const state = atlasSessionReducer(
      {
        ...initialAtlasSession,
        trajectory: {
          ...initialAtlasSession.trajectory,
          query: 'Astro',
          umapViewState: { target: [1, 2, 0], zoom: 3 },
        },
      },
      { type: 'openTrajectory', record },
    )
    expect(state.trajectory).toEqual({
      query: '',
      selectedId: 'HK1-Astro',
      umapViewState: { target: [1, 2, 0], zoom: 3 },
      umapResetVersion: 0,
    })
  })
})
