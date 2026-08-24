// 13 - Validate and publish the report package.
process BUILD_REPORT_PACKAGE {
    tag "$ligand_id"
    cpus 1; memory '512 MB'; time '5m'
    publishDir { "${params.outdir}/${ligand_id}" }, mode: params.publish_mode

    input:
    tuple val(ligand_id), path(params_json), path(report_files)

    output:
    tuple val(ligand_id), path('params.json'), path('package/')

    script:
    """
    python "${projectDir}/bin/build_report_package.py" --params params.json \
        --report-dir "$report_files" --out-dir package
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" build-package --ligand-id '${ligand_id}'
    """
}
