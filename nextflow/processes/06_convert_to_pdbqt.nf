// 06 - Prepare the ligand in PDBQT format.
process CONVERT_TO_PDBQT {
    tag "$ligand_id"
    cpus 1; memory '1 GB'; time '10m'

    input:
    tuple val(ligand_id), path(params_json), path(ligand_sdf), path(ligand_svg)

    output:
    tuple val(ligand_id), path('params.json'), path('ligand_input.sdf'),
          path('ligand_2d.svg'), path('ligand.pdbqt')

    script:
    """
    python "${projectDir}/bin/convert_ligand_to_pdbqt.py" --params params.json \
        --ligand-sdf ligand_input.sdf --out ligand.pdbqt
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" convert-pdbqt
    """
}
