// Validate the complete set of paths before dispatching analysis or rendering.
def parseCellTypeSelection(value) {
    if (value == null) return []
    def selected = value.toString().split(',').collect { part -> part.trim() }.findAll { part -> part }.unique()
    if (!selected) {
        error('--cell_types must contain at least one comma-separated cell-type slug')
    }
    return selected
}

// Let channel operators raise errors directly so Nextflow prints the useful message.
def selectCellTypes(entries, selected) {
    def slugs = entries.collect { entry -> entry[0].toString() }
    def duplicates = slugs.countBy { slug -> slug }.findAll { _slug, count -> count > 1 }.keySet()
    if (duplicates) {
        return [error: "Duplicate cell-type inputs: ${duplicates.sort().join(', ')}"]
    }
    if (!slugs) {
        return [error: 'No cell-type inputs are available']
    }
    def missing = selected - slugs
    if (missing) {
        return [error: (
            "Requested cell types are absent from the input: ${missing.sort().join(', ')}. " +
            "Available: ${slugs.sort().join(', ')}"
        )]
    }
    return [entries: selected ? entries.findAll { entry -> selected.contains(entry[0].toString()) } : entries]
}

def pairCellTypeInputs(entries, selected) {
    def splits = selectCellTypes(entries.findAll { entry -> entry[0] == 'split' }.collect { entry -> entry.drop(1) }, selected)
    def analyses = selectCellTypes(entries.findAll { entry -> entry[0] == 'analysis' }.collect { entry -> entry.drop(1) }, selected)
    if (splits.error) return splits
    if (analyses.error) return analyses
    def splitSlugs = splits.entries.collect { entry -> entry[0].toString() } as Set
    def analysisSlugs = analyses.entries.collect { entry -> entry[0].toString() } as Set
    if (splitSlugs != analysisSlugs) {
        return [error: (
            "Cell-type split/analysis inputs do not match. Missing analyses: ${(splitSlugs - analysisSlugs).sort()}; " +
            "missing splits: ${(analysisSlugs - splitSlugs).sort()}"
        )]
    }
    def bySlug = analyses.entries.collectEntries { entry -> [(entry[0].toString()): entry[1]] }
    return [entries: splits.entries.collect { entry -> [entry[0], entry[1], bySlug[entry[0].toString()]] }]
}
