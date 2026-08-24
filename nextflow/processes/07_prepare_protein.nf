// 07 - Repair and prepare the selected receptor.
process PREPARE_PROTEIN {
    tag "$ligand_id"
    cpus 2; memory '2 GB'; time '20m'
    errorStrategy { sleep(30 * task.attempt as long); 'retry' }
    maxRetries 3
    cache 'deep'

    input:
    tuple val(ligand_id), path(params_json), path(ligand_sdf),
          path(ligand_svg), path(ligand_pdbqt)

    output:
    tuple val(ligand_id), path('params.json'), path('ligand_input.sdf'),
          path('ligand_2d.svg'), path('ligand.pdbqt'), path('prepared_receptor.pdb'),
          path('receptor.pdbqt'), path('receptor.json'), path('protein_stats.json'),
          path('receptor_source.cif'), path('receptor_prep.log')

    script:
    """
    python "${projectDir}/bin/prepare_protein.py" --params params.json
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" prepare-protein
    """
}
