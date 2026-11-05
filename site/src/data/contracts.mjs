/** @typedef {{ gene: string; effect: number | null; padj: number | null; pvalue: number | null; se: number | null; baseMean: number | null; studySupport?: { significant: number; tested: number } | null }} Gene */
/** @typedef {{ id: string; model: string; study: string | null; covariate: string; contrast: string | null; design: string; samples: number; studies: Record<string, number>; records: Gene[]; recordsAsset?: string; recordBuckets?: Record<string, string>; downloadAsset?: string }} Model */
/** @typedef {{ id: string; label: string; tissue: string; analysisMode: string; warning: string | null; trajectoryAsset: string | null; referenceExpressionAsset?: string; previewAsset?: string; models: Model[] }} Cell */
/** @typedef {{ schemaVersion: number; releaseId?: string; source: string; cellTypes: Cell[] }} DifferentialRelease */
/** @typedef {{ schemaVersion: number; releaseId?: string; cellType: string; referenceBin: string; method: string; genes: Record<string, { cpm: number; studies: number; samples: number }> }} ReferenceExpression */
/** @typedef {{ id: string; x: number; y: number; cellType: string; tissue: string; gene: string; adjustedPValue: number; referenceExpression: number | null; referenceExpressionStudies?: number | null; referenceExpressionSamples?: number | null; meanCellsPerDonor: number | null; valueAtMaxAbsChange: number; cluster: string; trajectory: number[] }} TrajectoryRecord */
/** @typedef {{ id: string; x: number; y: number; cellType: string; tissue: string; gene: string; adjustedPValue: number; cluster: string }} TrajectoryPoint */
/** @typedef {{ release: string; warning?: string | null; description?: string; records: TrajectoryRecord[]; tissues: string[]; cellTypes: string[]; genes: string[]; geneBuckets?: Record<string, string> }} Atlas */
/** @typedef {{ schemaVersion: number; releaseId?: string; release: string; warning?: string | null; description?: string; records: TrajectoryPoint[]; tissues: string[]; cellTypes: string[]; genes: string[]; geneBuckets?: Record<string, string> }} TrajectoryMap */
/** @typedef {{ cellType: string; parent: string; study: string; samples: number; meanFraction: number; residualSd: number | null; binomialSamplingSd: number | null; covariates: string; warning: string | null }} Fraction */
/** @typedef {{ cellType: string; parent: string; study: string; covariate: string; term: string; contrast: string; estimate: number; low: number; high: number; samples: number }} FractionCoefficient */
/** @typedef {{ cellType: string; denominator: string; referenceBin: string; referenceMeanFraction: number | null; referenceAdjustedFraction: number | null; ageBin: string; group: string; estimate: number; low: number; high: number; samples: number }} FractionAgeEstimate */
/** @typedef {{ schemaVersion: number; records: Fraction[]; coefficients: FractionCoefficient[]; ageEstimates: FractionAgeEstimate[]; ageStatus: { cellType: string; study?: string; denominator?: string; status: string; reason?: string }[] }} FractionIndex */
/** @typedef {{ bin: string; effect: number | null }} Bin */
/** @typedef {{ referenceBin: string | null; design: string | null; genes: Record<string, Bin[]>; geneBuckets?: Record<string, string> }} Preview */
/** @typedef {{ id: string; label: string; description: string; genes: string[]; sourceFile: string }} PathwaySet */
/** @typedef {{ id: string; label: string }} GeneAnnotationTerm */
/** @typedef {{ families: GeneAnnotationTerm[]; proteinClasses: GeneAnnotationTerm[]; goTerms: { biologicalProcesses: GeneAnnotationTerm[]; molecularFunctions: GeneAnnotationTerm[]; cellularComponents: GeneAnnotationTerm[] } }} GeneMetadata */
/** @typedef {{ release: string; source: string; sourceUrl: string; sourceSha256: string; genes: Record<string, GeneMetadata> }} GeneMetadataIndex */
/** @typedef {{ id: string; label: string; path: string }} SingleCellType */
/** @typedef {{ schemaVersion: number; cellTypes: SingleCellType[] }} SingleCellIndex */
/** @typedef {{ gene: string; effect: number; padj: number | null }} ClusterMarker */
/** @typedef {{ cluster: string; image: string; fits: { study: string | null; samples: number; meanFraction: number; slopePpPerDecade: number | null }[] }} ClusterFractionDetail */
/** @typedef {{ schemaVersion: number; cellType: string; nCells: number; nClusters: number; integration: string; image: string; fractionImage?: string | null; fractionNote?: string | null; fractionHasStudyTrends?: boolean; fractionIncludesZeros?: boolean; fractionDetails?: ClusterFractionDetail[] | null; groups: { cluster: string; higher?: ClusterMarker[]; lower?: ClusterMarker[]; markers?: { gene: string; score: number; padj: number | null }[] }[] }} SingleCell */
/** @typedef {{ de: DifferentialRelease; trajectory: Atlas; trajectoryMap: TrajectoryMap; composition: FractionIndex; preview: Preview; pathways: { pathwaySets: PathwaySet[] }; metadata: GeneMetadataIndex; singleCell: SingleCell; singleCellIndex: SingleCellIndex; referenceExpression: ReferenceExpression }} Assets */

// Shared runtime contract for published JSON and browser loaders.
const obj = (x) => x !== null && typeof x === 'object' && !Array.isArray(x)
const str = (x) => typeof x === 'string' && x.length > 0
const nullableString = (x) => x === null || str(x)
const num = (x) => typeof x === 'number' && Number.isFinite(x)
const nullableNumber = (x) => x === null || num(x)
const count = (x) => Number.isInteger(x) && x >= 0
const probability = (x) => num(x) && x >= 0 && x <= 1
const arr = (x, predicate) => Array.isArray(x) && x.every(predicate)
const base = (x) => obj(x) && x.schemaVersion === 1
const gene = (x) =>
  obj(x) &&
  str(x.gene) &&
  nullableNumber(x.effect) &&
  (x.padj === null || probability(x.padj)) &&
  (x.pvalue === null || probability(x.pvalue)) &&
  nullableNumber(x.se) &&
  nullableNumber(x.baseMean) &&
  (x.studySupport == null || (obj(x.studySupport) && count(x.studySupport.significant) &&
    count(x.studySupport.tested) && x.studySupport.tested > 0 && x.studySupport.significant <= x.studySupport.tested))
const model = (x) =>
  obj(x) &&
  str(x.id) &&
  str(x.model) &&
  nullableString(x.study) &&
  str(x.covariate) &&
  nullableString(x.contrast) &&
  str(x.design) &&
  count(x.samples) &&
  obj(x.studies) &&
  Object.values(x.studies).every(count) &&
  (x.recordsAsset === undefined || str(x.recordsAsset)) &&
  (x.downloadAsset === undefined || str(x.downloadAsset)) &&
  (x.recordBuckets === undefined || (obj(x.recordBuckets) && Object.values(x.recordBuckets).every(str))) &&
  arr(x.records, gene)
const cell = (x) =>
  obj(x) &&
  str(x.id) &&
  str(x.label) &&
  str(x.tissue) &&
  str(x.analysisMode) &&
  nullableString(x.warning) &&
  nullableString(x.trajectoryAsset) &&
  (x.referenceExpressionAsset === undefined || str(x.referenceExpressionAsset)) &&
  (x.previewAsset === undefined || str(x.previewAsset)) &&
  arr(x.models, model) &&
  x.models.length > 0
const trajectory = (x) =>
  obj(x) &&
  str(x.id) &&
  num(x.x) &&
  num(x.y) &&
  str(x.cellType) &&
  str(x.tissue) &&
  str(x.gene) &&
  probability(x.adjustedPValue) &&
  nullableNumber(x.referenceExpression) &&
  (x.referenceExpressionStudies === undefined ||
    x.referenceExpressionStudies === null ||
    count(x.referenceExpressionStudies)) &&
  (x.referenceExpressionSamples === undefined ||
    x.referenceExpressionSamples === null ||
    count(x.referenceExpressionSamples)) &&
  nullableNumber(x.meanCellsPerDonor) &&
  num(x.valueAtMaxAbsChange) &&
  str(x.cluster) &&
  arr(x.trajectory, num) &&
  x.trajectory.length === 7
const fraction = (x) =>
  obj(x) &&
  str(x.cellType) &&
  str(x.parent) &&
  str(x.study) &&
  count(x.samples) &&
  probability(x.meanFraction) &&
  nullableNumber(x.residualSd) &&
  nullableNumber(x.binomialSamplingSd) &&
  typeof x.covariates === 'string' &&
  nullableString(x.warning)
const coefficient = (x) =>
  obj(x) &&
  ['cellType', 'parent', 'study', 'covariate', 'term', 'contrast'].every((k) => str(x[k])) &&
  ['estimate', 'low', 'high'].every((k) => num(x[k])) &&
  x.low <= x.high &&
  count(x.samples)
const age = (x) =>
  obj(x) &&
  ['cellType', 'group', 'denominator', 'ageBin', 'referenceBin'].every((k) => str(x[k])) &&
  ['estimate', 'low', 'high'].every((k) => num(x[k])) &&
  x.low <= x.high &&
  count(x.samples) &&
  nullableNumber(x.referenceMeanFraction) &&
  nullableNumber(x.referenceAdjustedFraction)
const term = (x) => obj(x) && str(x.id) && str(x.label)
const geneMetadata = (x) =>
  obj(x) &&
  arr(x.families, term) &&
  arr(x.proteinClasses, term) &&
  obj(x.goTerms) &&
  ['biologicalProcesses', 'molecularFunctions', 'cellularComponents'].every((k) => arr(x.goTerms[k], term))
const validators = {
  singleCellIndex: (x) => base(x) && arr(x.cellTypes, (c) => obj(c) && str(c.id) && str(c.label) && str(c.path)),
  singleCell: (x) => base(x) && str(x.cellType) && count(x.nCells) && count(x.nClusters) &&
    str(x.integration) && str(x.image) && x.image.startsWith('iVBORw0KGgo') &&
    (x.fractionImage === undefined || x.fractionImage === null || (str(x.fractionImage) && x.fractionImage.startsWith('iVBORw0KGgo'))) &&
    (x.fractionNote === undefined || x.fractionNote === null || typeof x.fractionNote === 'string') &&
    (x.fractionHasStudyTrends === undefined || typeof x.fractionHasStudyTrends === 'boolean') &&
    (x.fractionIncludesZeros === undefined || typeof x.fractionIncludesZeros === 'boolean') &&
    (x.fractionDetails == null || arr(x.fractionDetails, (d) => obj(d) && str(d.cluster) &&
      str(d.image) && d.image.startsWith('iVBORw0KGgo') && arr(d.fits, (fit) => obj(fit) &&
        nullableString(fit.study) && count(fit.samples) && probability(fit.meanFraction) &&
        nullableNumber(fit.slopePpPerDecade)))) &&
    arr(x.groups, (g) => obj(g) && str(g.cluster) && (
      (arr(g.higher, (m) => obj(m) && str(m.gene) && num(m.effect) && m.effect > 0 && (m.padj === null || probability(m.padj))) &&
       arr(g.lower, (m) => obj(m) && str(m.gene) && num(m.effect) && m.effect < 0 && (m.padj === null || probability(m.padj)))) ||
      arr(g.markers, (m) => obj(m) && str(m.gene) && num(m.score) && m.score > 0 && (m.padj === null || probability(m.padj)))
    )),
  de: (x) => base(x) && str(x.source) && arr(x.cellTypes, cell) && x.cellTypes.length > 0,
  referenceExpression: (x) =>
    base(x) &&
    str(x.cellType) &&
    str(x.referenceBin) &&
    str(x.method) &&
    obj(x.genes) &&
    Object.values(x.genes).every(
      (record) =>
        obj(record) &&
        num(record.cpm) &&
        record.cpm >= 0 &&
        count(record.studies) &&
        record.studies > 0 &&
        count(record.samples) &&
        record.samples >= record.studies,
    ),
  trajectory: (x) =>
    base(x) &&
    str(x.release) &&
    (x.warning === undefined || nullableString(x.warning)) &&
    (x.description === undefined || typeof x.description === 'string') &&
    (x.geneBuckets === undefined || (obj(x.geneBuckets) && Object.values(x.geneBuckets).every(str))) &&
    arr(x.records, trajectory) &&
    ['tissues', 'cellTypes', 'genes'].every((k) => arr(x[k], str)),
  trajectoryMap: (x) =>
    base(x) &&
    str(x.release) &&
    (x.warning === undefined || nullableString(x.warning)) &&
    (x.description === undefined || typeof x.description === 'string') &&
    (x.geneBuckets === undefined || (obj(x.geneBuckets) && Object.values(x.geneBuckets).every(str))) &&
    arr(x.records, (record) =>
      obj(record) && str(record.id) && num(record.x) && num(record.y) &&
      str(record.cellType) && str(record.tissue) && str(record.gene) &&
      probability(record.adjustedPValue) && str(record.cluster),
    ) &&
    ['tissues', 'cellTypes', 'genes'].every((k) => arr(x[k], str)),
  composition: (x) =>
    base(x) &&
    arr(x.records, fraction) &&
    arr(x.coefficients, coefficient) &&
    arr(x.ageEstimates, age) &&
    arr(
      x.ageStatus,
      (s) =>
        obj(s) &&
        str(s.cellType) &&
        str(s.status) &&
        (s.study === undefined || typeof s.study === 'string') &&
        (s.denominator === undefined || typeof s.denominator === 'string') &&
        (s.reason === undefined || typeof s.reason === 'string'),
    ),
  preview: (x) =>
    base(x) &&
    nullableString(x.referenceBin) &&
    nullableString(x.design) &&
    obj(x.genes) &&
    (x.geneBuckets === undefined || (obj(x.geneBuckets) && Object.values(x.geneBuckets).every(str))) &&
    Object.values(x.genes).every((bins) =>
      arr(bins, (b) => obj(b) && str(b.bin) && nullableNumber(b.effect)),
    ),
  pathways: (x) =>
    base(x) &&
    arr(
      x.pathwaySets,
      (p) =>
        obj(p) && str(p.id) && str(p.label) && str(p.description) && str(p.sourceFile) && arr(p.genes, str),
    ),
  metadata: (x) =>
    base(x) &&
    ['release', 'source', 'sourceUrl', 'sourceSha256'].every((k) => str(x[k])) &&
    obj(x.genes) &&
    Object.values(x.genes).every(geneMetadata),
}
/**
 * @template {keyof Assets} K
 * @param {K} kind
 * @param {unknown} value
 * @returns {value is Assets[K]}
 */
export function validAsset(kind, value) {
  return !!validators[kind]?.(value)
}
