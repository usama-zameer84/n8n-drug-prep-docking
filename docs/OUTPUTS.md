# Output reference

Every successful run creates one report folder and one `MD_Handoff` subfolder in Google Drive.

## Report folder

### `index.html`

The landing page gives the ligand and receptor identity, run ID, best Vina score, selected replica, QC flags, and a plain statement of what the run does and does not establish.

### `01_input_provenance.html`

Records the source Drive file ID and name, source checksum, input profile, original and canonical SMILES, RCSB entry, RCSB download URL, and receptor source SHA-256.

### `02_ligand_preparation.html`

Records the conformer method, random seed, force field, optimization status, Meeko version, charge method, atom mapping status, molecular descriptors, Rule-of-Five results, and the ligand's 2D depiction.

### `03_receptor_preparation.html`

Records available and selected chains, heterogen policy, detected waters and non-water components, missing residues, repaired atoms, pH, receptor atom counts, Meeko options, and the component used as the binding-site reference.

### `03a_receptor_selection.html`

Records automatic or supplied mode, ranked target findings, ChEMBL identifiers, UniProt accessions, assay evidence, maximum similarity, pChEMBL values, candidate PDB structures, resolution, author-chain mapping, preparation preflight outcomes, the selected receptor, and the reasons for selection.

### `04_docking_results.html`

Records the Vina version, box center and size, box source, exhaustiveness, requested modes, energy range, CPU setting, replica seeds, timeout, ranked poses, top scores, mean, standard deviation, range, and the score-interpretation warning.

### `05_distance_contacts.html`

Lists receptor residues within the configured distance cutoff of the docked ligand. Each row includes minimum distance, contact count, and closest geometric category. This is a distance screen, not a force-field interaction energy calculation.

### `06_methods_qc.html`

Collects software versions, reproducibility parameters, workflow QC flags, and the review rules for box placement, heterogens, metals, waters, protonation, missing residues, and replica behavior.

### `07_md_handoff_readiness.html`

States whether structural export passed atom-count checks, identifies the chemistry-authoritative ligand file, and lists the force-field, topology, solvation, equilibration, production, and convergence work still required.

### `08_raw_data.html`

Contains the complete report object in readable JSON. Use `run_summary.json` or files from `MD_Handoff` for programmatic pipelines.

### `09_visualization.html`

Loads the prepared receptor and every best pose in the matching receptor-profile registry. Controls show or hide individual ligands. The table records ligand name, run ID, color, Vina score, and box source.

The viewer needs network access to load 3Dmol.js. The SDF and PDB files remain the portable source data.

### Overlay artifacts

| File | Contents |
| --- | --- |
| `receptor_overlay_reference.pdb` | Exact prepared receptor used for the overlay group |
| `all_best_poses.sdf` | Concatenated best-pose SDF records |
| `ligand_overlay.json` | Receptor profile and per-pose metadata |
| `run_summary.json` | Compact run and Drive delivery summary |

The overlay key includes the RCSB source hash, selected chains, heterogen policy, and pH. Structures with a different preparation profile are placed in a different group.

## MD handoff folder

### Receptor files

- `receptor_source.cif`
- `prepared_receptor.pdb`
- `receptor.pdbqt`
- `receptor.json`
- `receptor_prep.log`
- `protein_stats.json`

### Ligand and pose files

- `ligand_input.sdf`
- `ligand.pdbqt`
- `ligand_2d.svg`
- `docked.pdbqt`
- `docked_poses.sdf`
- `best_pose_ligand.sdf`
- `best_pose_ligand.pdb`
- `complex_best_pose.pdb`
- `atom_mapping.json`

`best_pose_ligand.sdf` is the chemistry-authoritative ligand handoff. PDB files provide coordinates but do not preserve small-molecule bond orders as reliably as SDF.

### Results and provenance

- `docking.log`
- `docking_results.json`
- `interactions.json`
- `meeko_export.log`
- `provenance.json`
- `README_MD_HANDOFF.md`
- `manifest.json`
- `docked_replica_<n>.pdbqt`
- `vina_replica_<n>.log`

`manifest.json` records file sizes and SHA-256 checksums. `provenance.json` records ligand, receptor, preparation, box, docking, selection, and software information.

## What must be added before MD

The handoff deliberately stops before topology generation. A real MD setup still needs:

1. visual and chemical inspection of the selected pose;
2. a protein force field and compatible water model;
3. ligand atom typing, charges, and parameters;
4. reconciliation of atom names through `atom_mapping.json`;
5. decisions for termini, disulfides, cofactors, metals, retained waters, and protonation;
6. a periodic box, solvent, counterions, and target ionic strength;
7. minimization and staged NVT/NPT equilibration;
8. replicated production trajectories and convergence analysis.
