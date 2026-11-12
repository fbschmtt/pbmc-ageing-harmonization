import { expect, it } from 'vitest'
import { navigation, navigationItems } from './navigation'

it('groups expression and age-analysis routes while keeping pathways as one page', () => {
  const ids = navigation.map((item) => item.id)
  const labels = navigationItems.filter((item) => 'pages' in item).map((item) => item.label)
  expect(ids).toContain('crossType')
  expect(navigationItems.some((item) => 'pages' in item && item.pages.some((page) => page.id === 'crossType'))).toBe(true)
  expect(ids).toContain('pathways')
  expect(ids).toContain('geneTypeMatrix')
  expect(ids.filter((id) => id === 'geneTypeMatrix')).toHaveLength(1)
  expect(labels).toContain('Differential expression')
  expect(labels).not.toContain('Pathways')
  expect(navigationItems.find((item) => item.id === 'pathways')?.label).toBe('Pathways')
  expect(navigationItems.some((item) => 'pages' in item && item.pages.some((page) => page.id === 'geneTypeMatrix'))).toBe(false)
  expect(labels).toContain('More age analysis')
  expect(navigation.find((item) => item.id === 'trajectories')?.label).toBe('Age-bin trajectories')
  expect(navigation.find((item) => item.id === 'genes')?.label).toBe('By gene')
})

it('includes the clock placeholder and about page without a download tab', () => {
  const ids = navigation.map((item) => item.id)
  expect(ids).toContain('clock')
  expect(ids).toContain('about')
  expect(ids).not.toContain('download')
  expect(navigation.find((item) => item.id === 'clock')?.label).toBe('Transcriptomic clocks')
})
