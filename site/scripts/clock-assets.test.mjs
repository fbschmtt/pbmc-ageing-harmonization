import { describe, expect, it } from 'vitest'
import { clockModelAssetPath } from './clock-assets.mjs'

describe('clock model asset paths', () => {
  it('maps cross-type method IDs to safe, deterministic paths', () => {
    expect(clockModelAssetPath('cross-cell-type-mean', 'cross_type_mean:elastic_net'))
      .toBe('clock-models/cross-cell-type-mean/cross_type_mean_elastic_net.json')
  })
})
