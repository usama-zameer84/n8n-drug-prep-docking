// 08 - Run the configured AutoDock Vina replicas.
process RUN_DOCKING {
    tag "$ligand_id"
    cpus params.cpu as int; memory '2 GB'; time "${params.timeout_seconds + 120}s"
    errorStrategy 'terminate'

    input:
    tuple val(ligand_id), path(params_json), path(ligand_sdf), path(ligand_svg),
          path(ligand_pdbqt), path(prepared_receptor), path(receptor_pdbqt),
          path(receptor_json), path(protein_stats), path(receptor_source), path(receptor_prep_log)

    output:
    tuple val(ligand_id), path('params.json'), path('ligand_input.sdf'), path('ligand_2d.svg'),
          path('ligand.pdbqt'), path('prepared_receptor.pdb'), path('receptor.pdbqt'),
          path('receptor.json'), path('protein_stats.json'), path('receptor_source.cif'),
          path('receptor_prep.log'), path('docked.pdbqt'), path('docking.log'),
          path('docked_replica_*.pdbqt'), path('vina_replica_*.log')

    script:
    """
    python "${projectDir}/bin/run_docking.py" --params params.json
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" run-docking
    """
}
