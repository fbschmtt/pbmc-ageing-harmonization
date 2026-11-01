import { expect, it } from 'vitest'
import { buildStudyColors } from './studyColors'

it('assigns colors consistently regardless of release study ordering or duplicates', () => {
  const first = buildStudyColors(['study_b', 'study_a', 'study_c'])
  const second = buildStudyColors(['study_c', 'study_a', 'study_b', 'study_a'])
  expect([...first]).toEqual([...second])
  expect(new Set(first.values()).size).toBe(3)
})
