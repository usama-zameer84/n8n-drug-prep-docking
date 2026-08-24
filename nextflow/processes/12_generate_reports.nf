// 12 - Build report data and render the report templates.
process GENERATE_REPORTS {
    tag "$ligand_id"
    cpus 1; memory '1 GB'; time '5m'

    input:
    tuple val(ligand_id), path(params_json), path(md_handoff_dir)

    output:
    tuple val(ligand_id), path('params.json'), path('report_files/')

    script:
    """
    python "${projectDir}/bin/generate_reports.py" --params params.json --out-dir report_files \
        --templates "${projectDir}/templates" --handoff-dir "$md_handoff_dir"
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" generate-reports --input-dir "$md_handoff_dir"
    """
}
