// 04 - Generate a three-dimensional ligand conformer.
process GENERATE_3D_STRUCTURE {
    tag "$ligand_id"
    cpus 1; memory '1 GB'; time '10m'

    input:
    tuple val(ligand_id), path(params_json)

    output:
    tuple val(ligand_id), path('params.json'), path('ligand_input.sdf'), path('ligand_2d.svg')

    script:
    """
    python "${projectDir}/bin/generate_3d_structure.py" --params params.json \
        --ligand-sdf ligand_input.sdf --ligand-svg ligand_2d.svg
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" generate-3d
    """
}
