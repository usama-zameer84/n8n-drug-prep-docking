// 03 - Select and preflight a receptor structure.
process AUTO_SELECT_RECEPTOR {
    tag "$ligand_id"
    cpus 1; memory '1 GB'; time '20m'
    errorStrategy { sleep(30 * task.attempt as long); 'retry' }
    maxRetries 3
    cache 'deep'

    input:
    tuple val(ligand_id), path(params_json)

    output:
    tuple val(ligand_id), path('params.json')

    script:
    """
    python "${projectDir}/bin/select_receptor.py" --params params.json
    """

    stub:
    """
    python "${projectDir}/bin/stub_stage.py" select-receptor
    """
}
