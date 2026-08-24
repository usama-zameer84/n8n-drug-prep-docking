// 10 - Find heavy-atom contacts for the selected pose.
process ANALYZE_INTERACTIONS {
    tag "$ligand_id"
    cpus 1; memory '1 GB'; time '5m'

    input:
    tuple val(ligand_id), path(params_json), path(ligand_sdf), path(ligand_svg),
          path(ligand_pdbqt), path(prepared_receptor), path(receptor_pdbqt),
          path(receptor_json), path(protein_stats), path(receptor_source), path(receptor_prep_log),
          path(docked_pdbqt), path(docking_log), path(docking_results),
          path(replica_poses), path(replica_logs)

    output:
    tuple val(ligand_id), path('params.json'), path('ligand_input.sdf'), path('ligand_2d.svg'),
          path('ligand.pdbqt'), path('prepared_receptor.pdb'), path('receptor.pdbqt'),
          path('receptor.json'), path('protein_stats.json'), path('receptor_source.cif'),
          path('receptor_prep.log'), path('docked.pdbqt'), path('docking.log'),
          path('docking_results.json'), path('interactions.json'), path(replica_poses), path(replica_logs)

    script:
    """
    python "${projectDir}/bin/analyze_interactions.py" --params params.json
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" analyze-interactions
    """
}
