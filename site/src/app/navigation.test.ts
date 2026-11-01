import { expect, it } from 'vitest'
import { navigation } from './navigation'

it('keeps the shared pathway selector and age-bin matrix as distinct routes', () => {
  const ids = navigation.map((item) => item.id)
  expect(ids).toContain('pathways')
  expect(ids).toContain('geneTypeMatrix')
  expect(ids.filter((id) => id === 'geneTypeMatrix')).toHaveLength(1)
  expect(navigation.find((item) => item.id === 'trajectories')?.label).toBe('Cross-type trajectory map')
})

it('places datasets and modelling immediately before download', () => {
  expect(navigation.slice(-3).map((item) => item.id)).toEqual(['datasets', 'modelling', 'download'])
})
