export type SelectableModel = { id: string; model: string; covariate: string }
export type ModelGroup<T> = { label: string; models: T[] }

export function groupModels<T extends SelectableModel>(models: T[]): ModelGroup<T>[] {
  const groups = [
    {
      label: 'Combined',
      models: models.filter((model) => model.model !== 'per_study' && model.model !== 'developmental'),
    },
    { label: 'Developmental', models: models.filter((model) => model.model === 'developmental') },
    { label: 'Individual study fits', models: models.filter((model) => model.model === 'per_study') },
  ]
  return groups.filter((group) => group.models.length > 0)
}

export function defaultModelId(models: SelectableModel[]): string {
  const nonStudy = models.filter((model) => model.model !== 'per_study')
  return (
    (nonStudy.find((model) => model.model === 'merged' && model.covariate === 'age') || nonStudy[0])?.id ?? ''
  )
}

export function modelLabel(
  model: SelectableModel & { study: string | null; contrast: string | null },
): string {
  const display: Record<string, string> = {
    age: 'Age',
    sex: 'Sex',
    bmi: 'BMI',
    cmv: 'CMV',
    smoking_status: 'Smoking status',
    log10_total_counts: 'Log₁₀ total counts',
    adult_age: 'Adult age',
    juvenile_age: 'Juvenile age',
    age_25_vs_5: 'Age 25 vs 5',
  }
  const raw = model.covariate.replaceAll('_', ' ')
  const covariate = display[model.covariate] || raw.charAt(0).toUpperCase() + raw.slice(1)
  const group =
    model.model === 'per_study'
      ? model.study || 'Individual study'
      : model.model === 'developmental'
        ? 'Developmental'
        : 'Combined'
  const contrast = model.covariate === 'age' && model.contrast === 'per year'
    ? 'per decade'
    : model.contrast
  return `${group} · ${covariate}${contrast ? ` (${contrast})` : ''}`
}
