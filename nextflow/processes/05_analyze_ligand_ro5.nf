// 05 - Calculate ligand descriptors and Rule-of-Five status.
process ANALYZE_LIGAND_RO5 {
    tag "$ligand_id"
    cpus 1; memory '512 MB'; time '5m'

    input:
    tuple val(ligand_id), path(params_json), path(ligand_sdf), path(ligand_svg)

    output:
    tuple val(ligand_id), path('params.json'), path('ligand_input.sdf'), path('ligand_2d.svg')

    script:
    """
    python "${projectDir}/bin/analyze_ligand_ro5.py" --params params.json \
        --ligand-sdf ligand_input.sdf
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" analyze-ro5
    """
}
