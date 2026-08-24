// 11 - Build the structural MD handoff bundle.
process BUILD_MD_HANDOFF {
    tag "$ligand_id"
    cpus 1; memory '1 GB'; time '10m'

    input:
    tuple val(ligand_id), path(params_json), path(ligand_sdf), path(ligand_svg),
          path(ligand_pdbqt), path(prepared_receptor), path(receptor_pdbqt),
          path(receptor_json), path(protein_stats), path(receptor_source), path(receptor_prep_log),
          path(docked_pdbqt), path(docking_log), path(docking_results), path(interactions),
          path(replica_poses), path(replica_logs)

    output:
    tuple val(ligand_id), path('params.json'), path('MD_Handoff/')

    script:
    """
    python "${projectDir}/bin/build_md_handoff.py" --params params.json --out-dir MD_Handoff
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" build-md-handoff
    """
}
