import type { TrajectoryRecord } from '../data/atlas'
import type { PathwaySet } from '../data/pathways'
import type { ExplorerState } from '../features/explorer/Explorer'
import type { ProfileState, ProfileTab } from '../features/profiles/ProfilesPage'

export type AtlasSessionState = {
  geneProfile: ProfileState
  cellTypeProfile: ProfileState
  trajectory: ExplorerState
  selectedGene: string | null
  selectedCellType: string | null
  selectedDEModelId: string | null
  compositionDenominator: string
  highlightGenes: string
  activePathwayId: string | null
  geneSetTouched: boolean
}

export const initialAtlasSession: AtlasSessionState = {
  geneProfile: { gene: null, cellType: null, selectedId: null },
  cellTypeProfile: { gene: null, cellType: null, selectedId: null },
  trajectory: {
    query: '',
    selectedId: null,
    umapViewState: null,
    umapResetVersion: 0,
  },
  selectedGene: null,
  selectedCellType: null,
  selectedDEModelId: null,
  compositionDenominator: 'l1_parent',
  highlightGenes: '',
  activePathwayId: null,
  geneSetTouched: false,
}

export type AtlasSessionAction =
  | { type: 'initializePathway'; pathways: PathwaySet[] }
  | { type: 'selectGene'; gene: string | null; cellType?: string }
  | { type: 'selectCellType'; cellType: string }
  | { type: 'selectDEModel'; modelId: string }
  | { type: 'selectCompositionDenominator'; denominator: string }
  | { type: 'editGeneSet'; value: string }
  | { type: 'selectPathway'; pathway: PathwaySet }
  | { type: 'setProfile'; tab: ProfileTab; state: ProfileState }
  | { type: 'openProfile'; tab: ProfileTab; record: TrajectoryRecord; selectPair: boolean }
  | { type: 'setTrajectory'; state: ExplorerState }
  | { type: 'openTrajectory'; record: TrajectoryRecord }

export function atlasSessionReducer(state: AtlasSessionState, action: AtlasSessionAction): AtlasSessionState {
  switch (action.type) {
    case 'initializePathway': {
      if (state.geneSetTouched || action.pathways.length === 0) return state
      const pathway =
        action.pathways.find((item) => item.id === 'REACTOME_CELLULAR_SENESCENCE') ?? action.pathways[0]
      return {
        ...state,
        highlightGenes: pathway.genes.join('\n'),
        activePathwayId: pathway.id,
        geneSetTouched: true,
      }
    }
    case 'selectGene': {
      const changed =
        action.gene !== state.selectedGene ||
        (action.cellType != null && action.cellType !== state.selectedCellType)
      if (!changed) return state
      return {
        ...state,
        selectedGene: action.gene,
        selectedCellType: action.cellType ?? state.selectedCellType,
        geneProfile: { ...state.geneProfile, gene: action.gene, selectedId: null },
        cellTypeProfile: {
          ...state.cellTypeProfile,
          cellType: action.cellType ?? state.cellTypeProfile.cellType,
          selectedId: null,
        },
        trajectory: { ...state.trajectory, query: '', selectedId: null },
      }
    }
    case 'selectCellType':
      if (action.cellType === state.selectedCellType) return state
      return {
        ...state,
        selectedCellType: action.cellType,
        geneProfile: { ...state.geneProfile, selectedId: null },
        cellTypeProfile: { ...state.cellTypeProfile, cellType: action.cellType, selectedId: null },
        trajectory: { ...state.trajectory, selectedId: null },
      }
    case 'selectDEModel':
      return state.selectedDEModelId === action.modelId
        ? state
        : { ...state, selectedDEModelId: action.modelId }
    case 'selectCompositionDenominator':
      return state.compositionDenominator === action.denominator
        ? state
        : { ...state, compositionDenominator: action.denominator }
    case 'editGeneSet':
      return { ...state, highlightGenes: action.value, activePathwayId: null, geneSetTouched: true }
    case 'selectPathway':
      return {
        ...state,
        highlightGenes: action.pathway.genes.join('\n'),
        activePathwayId: action.pathway.id,
        geneSetTouched: true,
      }
    case 'setProfile':
      return action.tab === 'gene'
        ? { ...state, geneProfile: action.state }
        : { ...state, cellTypeProfile: action.state }
    case 'openProfile': {
      const profile = {
        gene: action.record.gene,
        cellType: action.record.cellType,
        selectedId: action.selectPair ? action.record.id : null,
      }
      return action.tab === 'gene'
        ? { ...state, geneProfile: profile }
        : { ...state, cellTypeProfile: profile }
    }
    case 'setTrajectory':
      return { ...state, trajectory: action.state }
    case 'openTrajectory':
      return {
        ...state,
        trajectory: {
          query: '',
          selectedId: action.record.id,
          umapViewState: null,
          umapResetVersion: 0,
        },
      }
  }
}
