nextflow.enable.dsl = 2

include { parseCellTypeSelection; selectCellTypes; pairCellTypeInputs } from '../../modules/cell_type_selection'

workflow {
    def inputs = [['cd14-monocyte', 'monocyte.h5ad'], ['naive-cd4-t-cell', 't-cell.h5ad']]
    assert parseCellTypeSelection(null) == []
    assert parseCellTypeSelection(' cd14-monocyte,cd14-monocyte, naive-cd4-t-cell ') == inputs.collect { entry -> entry[0] }
    assert selectCellTypes(inputs, []).entries == inputs
    def pairs = [
        ['split', 'cd14-monocyte', 'monocyte.h5ad'],
        ['analysis', 'cd14-monocyte', 'monocyte-analysis'],
        ['split', 'naive-cd4-t-cell', 't-cell.h5ad'],
    ]
    assert pairCellTypeInputs(pairs, ['cd14-monocyte']).entries == [
        ['cd14-monocyte', 'monocyte.h5ad', 'monocyte-analysis'],
    ]
    // Exercise the same collect/flatMap boundary used by the production workflow.
    channel.fromList(inputs).collect(flat: false)
        .flatMap { entries -> selectCellTypes(entries, ['cd14-monocyte']).entries }
        .collect(flat: false)
        .view { selected ->
            assert selected == [inputs[0]]
            'Cell-type selection and input pairing checks passed'
        }
}
