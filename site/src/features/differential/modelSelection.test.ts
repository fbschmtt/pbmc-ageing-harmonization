import { describe, expect, it } from 'vitest'
import { defaultModelId, groupModels, modelLabel } from './modelSelection'

const combinedAge = {
  id: 'combined-age',
  model: 'merged',
  covariate: 'age',
  study: null,
  contrast: 'per year',
}
const combinedSex = {
  id: 'combined-sex',
  model: 'combined',
  covariate: 'sex',
  study: null,
  contrast: 'male vs female',
}
const developmental = {
  id: 'developmental-age',
  model: 'developmental',
  covariate: 'adult_age',
  study: null,
  contrast: null,
}
const studyAge = {
  id: 'study-age',
  model: 'per_study',
  covariate: 'age',
  study: 'study_a',
  contrast: 'per year',
}

describe('differential expression model selection', () => {
  it('groups combined and developmental fits before individual study fits', () => {
    const models = [studyAge, combinedSex, developmental, combinedAge]
    expect(groupModels(models)).toEqual([
      { label: 'Combined', models: [combinedSex, combinedAge] },
      { label: 'Developmental', models: [developmental] },
      { label: 'Individual study fits', models: [studyAge] },
    ])
    expect(defaultModelId(models)).toBe('combined-age')
  })
  it('keeps study fits available if no combined model is available', () => {
    expect(defaultModelId([studyAge])).toBe('')
    expect(groupModels([studyAge])[0].models).toEqual([studyAge])
  })
  it('uses consistent capitalization without changing study IDs or contrasts', () => {
    expect(modelLabel(combinedAge)).toBe('Combined · Age (per year)')
    expect(modelLabel(developmental)).toBe('Developmental · Adult age')
    expect(modelLabel(studyAge)).toBe('study_a · Age (per year)')
  })
})
