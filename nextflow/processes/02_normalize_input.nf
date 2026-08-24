// 02 - Validate the ligand input and write params.json.
process NORMALIZE_INPUT {
    tag "$ligand_id"
    cpus 1; memory '512 MB'; time '5m'

    input:
    tuple val(ligand_id), val(smiles)

    output:
    tuple val(ligand_id), path('params.json')

    script:
    def inputJson = groovy.json.JsonOutput.toJson([
        ligand_id: ligand_id,
        smiles: smiles,
        receptor_selection_mode: params.receptor_selection_mode,
        pdb_id: params.pdb_id,
        chain_ids: params.chain_ids,
        center_x: params.center_x, center_y: params.center_y, center_z: params.center_z,
        size_x: params.size_x, size_y: params.size_y, size_z: params.size_z,
        ph: params.ph, exhaustiveness: params.exhaustiveness,
        num_modes: params.num_modes, energy_range: params.energy_range,
        seed: params.seed, cpu: params.cpu, replicas: params.replicas,
        timeout_seconds: params.timeout_seconds, cutoff: params.cutoff,
        target_organism: params.target_organism,
        target_similarity_threshold: params.target_similarity_threshold,
        target_candidate_limit: params.target_candidate_limit,
        heterogen_policy: params.heterogen_policy,
        add_missing_residues: params.add_missing_residues,
        allow_multicomponent: params.allow_multicomponent,
    ])
    def inputBase64 = inputJson.bytes.encodeBase64().toString()
    """
    python "${projectDir}/bin/normalize_input.py" --input-base64 '${inputBase64}' --out params.json
    """

    stub:
    def stubJson = groovy.json.JsonOutput.toJson([ligand_id: ligand_id, smiles: smiles])
    def stubBase64 = stubJson.bytes.encodeBase64().toString()
    """
    python "${projectDir}/bin/stub_stage.py" normalize --input-base64 '${stubBase64}'
    """
}
