// Ligand preparation, receptor selection, docking, reports, and structural handoff.

include { NORMALIZE_INPUT }       from './processes/02_normalize_input.nf'
include { AUTO_SELECT_RECEPTOR }  from './processes/03_auto_select_receptor.nf'
include { GENERATE_3D_STRUCTURE } from './processes/04_generate_3d_structure.nf'
include { ANALYZE_LIGAND_RO5 }    from './processes/05_analyze_ligand_ro5.nf'
include { CONVERT_TO_PDBQT }      from './processes/06_convert_to_pdbqt.nf'
include { PREPARE_PROTEIN }       from './processes/07_prepare_protein.nf'
include { RUN_DOCKING }           from './processes/08_run_docking.nf'
include { PARSE_RESULTS }         from './processes/09_parse_results.nf'
include { ANALYZE_INTERACTIONS }  from './processes/10_analyze_interactions.nf'
include { BUILD_MD_HANDOFF }      from './processes/11_build_md_handoff.nf'
include { GENERATE_REPORTS }      from './processes/12_generate_reports.nf'
include { BUILD_REPORT_PACKAGE }  from './processes/13_build_report_package.nf'

workflow {
    def inputFile = file(params.input, checkIfExists: true)
    def ligandRecords = inputFile.readLines()
        .collect { line -> line.trim() }
        .findAll { line -> line && !line.startsWith('#') }
        .collect { line ->
            def parts  = line.split(/\s+/)
            def smiles = parts[0]
            def name   = parts.length > 1 ? parts[1..-1].join('_') : 'ligand'
            def safe   = name.replaceAll(/[^A-Za-z0-9_.-]+/, '_').replaceAll(/(^_+|_+$)/, '') ?: 'ligand'
            tuple(safe, smiles)
        }
    if (ligandRecords.size() != 1) {
        error "Input must contain exactly one ligand record: ${inputFile}"
    }
    ch_ligands = channel.of(ligandRecords[0])

    // 02 normalize: validate SMILES + params, emit params.json
    ch_norm = NORMALIZE_INPUT(ch_ligands)

    // 03 auto-select receptor (ChEMBL/RCSB/UniProt + PDBFixer/Meeko preflight)
    ch_rec  = AUTO_SELECT_RECEPTOR(ch_norm)

    // 04 generate 3D structure (rdkit)
    ch_3d   = GENERATE_3D_STRUCTURE(ch_rec)

    // 05 Lipinski Rule of 5
    ch_ro5  = ANALYZE_LIGAND_RO5(ch_3d)

    // 06 convert ligand to PDBQT (meeko mk_prepare_ligand)
    ch_pdbqt = CONVERT_TO_PDBQT(ch_ro5)

    // 07 prepare protein (PDBFixer + mk_prepare_receptor)
    ch_prot = PREPARE_PROTEIN(ch_pdbqt)

    // 08 run docking (AutoDock Vina)
    ch_dock = RUN_DOCKING(ch_prot)

    // 09 parse Vina results
    ch_parse = PARSE_RESULTS(ch_dock)

    // 10 analyze interactions
    ch_int  = ANALYZE_INTERACTIONS(ch_parse)

    // 11 build MD handoff bundle (structural-only)
    ch_md   = BUILD_MD_HANDOFF(ch_int)

    // 12 generate HTML reports
    ch_rep  = GENERATE_REPORTS(ch_md)

    // 13 build per-ligand report package + publish to results/<ligand_id>/
    ch_pkg  = BUILD_REPORT_PACKAGE(ch_rep)

}
