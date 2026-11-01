type Rgb = [number, number, number]
type Rgba = [number, number, number, number]

function mixHex(foreground: string, background: string, foregroundWeight: number) {
  const foregroundRgb = rgb(foreground)
  const backgroundRgb = rgb(background)
  return `#${foregroundRgb
    .map((channel, index) =>
      Math.round(channel * foregroundWeight + backgroundRgb[index] * (1 - foregroundWeight))
        .toString(16)
        .padStart(2, '0'),
    )
    .join('')}`
}

const palette = {
  canvas: '#f8f7f2',
  brain: '#256b6e',
  blood: '#d98c63',
} as const

export const visualColors = {
  ...palette,
  subduedBrain: mixHex(palette.brain, palette.canvas, 0.1),
  subduedBlood: mixHex(palette.blood, palette.canvas, 0.1),
  highlightedBrain: '#164f52',
  highlightedBlood: '#ad5f3e',
  geneSetOutline: '#164f52',
  selected: '#b5453c',
  positiveChange: '#bf634d',
  negativeChange: '#27786f',
  trajectoryRise: '#27786f',
  trajectoryFall: '#bd6e54',
  missingResult: '#ecefe9',
} as const

function rgb(hex: string): Rgb {
  return [
    Number.parseInt(hex.slice(1, 3), 16),
    Number.parseInt(hex.slice(3, 5), 16),
    Number.parseInt(hex.slice(5, 7), 16),
  ]
}

function rgba(hex: string, alpha = 255): Rgba {
  return [...rgb(hex), alpha]
}

export const visualRgb = {
  positiveChange: rgb(visualColors.positiveChange),
  negativeChange: rgb(visualColors.negativeChange),
} as const

export const deckColors = {
  subduedBrain: rgba(visualColors.subduedBrain),
  subduedBlood: rgba(visualColors.subduedBlood),
  brain: rgba(visualColors.brain, 184),
  blood: rgba(visualColors.blood, 184),
  highlightedBrain: rgba(visualColors.highlightedBrain),
  highlightedBlood: rgba(visualColors.highlightedBlood),
  selected: rgba(visualColors.selected),
} as const

export function installVisualEncoding(element: HTMLElement) {
  Object.entries(visualColors).forEach(([name, value]) => {
    const cssName = name.replace(/[A-Z]/g, (letter) => `-${letter.toLowerCase()}`)
    element.style.setProperty(`--atlas-${cssName}`, value)
  })
}
