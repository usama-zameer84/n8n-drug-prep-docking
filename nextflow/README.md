# Nextflow workflow

The Nextflow workflow prepares and docks one ligand per run. It accepts a SMILES file, selects
or validates a receptor, runs AutoDock Vina, publishes an HTML report package, and writes a
structural MD handoff.

Docking scores rank predicted poses. They do not measure experimental affinity or establish
binding. The MD handoff is not topology-ready.

## Directory contents

| Path | Contents |
| --- | --- |
| `main.nf` | Workflow entry point and single-ligand input parser |
| `nextflow.config` | Defaults, process settings, and execution profiles |
| `processes/` | Processes 02 through 13; orchestration and resource declarations only |
| `bin/` | Scientific, reporting, packaging, and stub-run Python helpers |
| `templates/` | HTML report templates and the offline structure viewer |
| `examples/input.smi` | Default caffeine input |
| `examples/input.test.smi` | Input used by the test profile and CI stub run |

Process definitions call helpers from `bin/`; they do not contain embedded Python programs.
This keeps workflow orchestration separate from the code that transforms scientific data and
builds report artifacts.

## Install

From the repository root:

```bash
conda env create -f environment.yml -n drugprep
conda activate drugprep

which nextflow vina
nextflow -version
vina --version
python -c "import rdkit, openmm, pdbfixer, meeko, gemmi, vina; print('scientific stack OK')"
```

The `local` profile uses the active environment. The `standard` and `conda` profiles ask
Nextflow to provision [`../environment.yml`](../environment.yml).

## Input file

The input is a text file containing exactly one non-comment ligand record:

```text
SMILES optional_ligand_name
```

Example:

```text
# caffeine
CN1C=NC2=C1C(=O)N(C(=O)N2C)C caffeine
```

Rules:

- blank lines and lines beginning with `#` are ignored
- the first whitespace-separated value is the SMILES
- remaining values form the ligand name, joined with underscores
- unsafe filename characters in the name are replaced with underscores
- no name produces the ID `ligand`
- zero records or more than one record stops the workflow
- the SMILES may contain at most 1000 characters
- dot-separated components require `--allow_multicomponent true`

`replicas` runs Vina more than once for this one ligand. It does not enable multiple ligand
records.

## Run examples

Use the active conda environment:

```bash
nextflow run nextflow/main.nf -profile local \
  --input nextflow/examples/input.smi
```

Let Nextflow create or reuse the environment:

```bash
nextflow run nextflow/main.nf -profile standard \
  --input nextflow/examples/input.smi
```

Use automatic receptor selection with stronger Vina search settings:

```bash
nextflow run nextflow/main.nf -profile local \
  --input nextflow/examples/input.smi \
  --target_organism "Homo sapiens" \
  --target_similarity_threshold 75 \
  --exhaustiveness 32 \
  --replicas 3
```

Use a supplied receptor and selected chains:

```bash
nextflow run nextflow/main.nf -profile local \
  --input nextflow/examples/input.smi \
  --receptor_selection_mode provided \
  --pdb_id 9H37 \
  --chain_ids A
```

Use a supplied docking box:

```bash
nextflow run nextflow/main.nf -profile local \
  --input nextflow/examples/input.smi \
  --center_x 12.5 --center_y 28.0 --center_z 40.0 \
  --size_x 20 --size_y 20 --size_z 20
```

Provide all three center coordinates or omit all three.

## Parameters

### Workflow parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `input` | `nextflow/examples/input.smi` | Single-ligand SMILES file |
| `outdir` | `<launch directory>/results` | Published result directory |
| `publish_mode` | `copy` | Nextflow `publishDir` mode |

### Receptor parameters

| Parameter | Default | Accepted value | Meaning |
| --- | --- | --- | --- |
| `receptor_selection_mode` | `auto` | `auto` or `provided` | Select a receptor or use `pdb_id` |
| `pdb_id` | empty | four-character PDB ID | Required for provided mode |
| `chain_ids` | empty | comma-separated chain IDs | Chains retained in provided mode |
| `target_organism` | `Homo sapiens` | non-empty name under 100 characters | Organism filter for activity evidence |
| `target_similarity_threshold` | `70` | 40 to 100 | Minimum ChEMBL similarity percentage |
| `target_candidate_limit` | `5` | integer from 1 to 10 | Target candidates retained for ranking |
| `heterogen_policy` | `remove_all` | `remove_all` or `keep_water` | Receptor heterogen handling |
| `add_missing_residues` | `false` | boolean | Allow PDBFixer to build missing residues |
| `ph` | `7.4` | 4 to 10 | Hydrogen-addition pH |

### Docking and analysis parameters

| Parameter | Default | Accepted value | Meaning |
| --- | --- | --- | --- |
| `center_x`, `center_y`, `center_z` | inferred | all three numbers from -10000 to 10000 Å | Docking-box center |
| `size_x`, `size_y`, `size_z` | `20` | 8 to 30 Å each; volume at most 27000 Å³ | Docking-box dimensions |
| `exhaustiveness` | `8` | integer from 8 to 64 | Vina search effort |
| `num_modes` | `9` | integer from 1 to 20 | Maximum poses retained per replica |
| `energy_range` | `3` | 1 to 10 kcal/mol | Vina score window |
| `seed` | `20260824` | integer from 1 to 2147483647 | Base random seed |
| `cpu` | `1` | integer from 1 to 8 | Vina CPU threads |
| `replicas` | `1` | integer from 1 to 3 | Repeated Vina runs for the same ligand |
| `timeout_seconds` | `900` | integer from 60 to 1800 | Timeout for each replica |
| `cutoff` | `4.5` | 2.5 to 6 Å | Heavy-atom contact cutoff |
| `allow_multicomponent` | `false` | boolean | Permit dot-separated components in the SMILES |

Without an explicit center, the workflow first uses the centroid of the most similar eligible
co-crystallized component. It falls back to the prepared-receptor centroid when necessary and
records that fallback in the QC flags.

## Execution profiles

| Profile | Use |
| --- | --- |
| `standard` | Local executor with a Nextflow-managed conda environment |
| `local` | Local executor using the already active environment |
| `conda` | Explicit local executor and Nextflow-managed conda environment |
| `slurm` | Shared-filesystem SLURM template; adjust queue and site resources |
| `awsbatch` | Configuration template; add a container, S3 paths, region, and real queue |
| `test` | Reduced settings used with another profile for smoke tests |

Profiles can be combined. For example, `-profile local,test` uses the active environment and
the smoke-test settings.

## Outputs

The default output tree is:

```text
results/<ligand_id>/
├── params.json
└── package/
    ├── index.html
    ├── 01_input_provenance.html
    ├── 02_ligand_preparation.html
    ├── 03_receptor_preparation.html
    ├── 03a_receptor_selection.html
    ├── 04_docking_results.html
    ├── 05_distance_contacts.html
    ├── 06_methods_qc.html
    ├── 07_md_handoff_readiness.html
    ├── 08_raw_data.html
    ├── 09_visualization.html
    ├── run_summary.json
    ├── ligand_overlay.json
    ├── all_best_poses.sdf
    ├── receptor_overlay_reference.pdb
    └── MD_Handoff/
```

`params.json` is the final run state passed between processes. It contains normalized input,
receptor-selection evidence, preparation metadata, docking results, contacts, QC flags,
handoff metadata, and the report file list.

### HTML and summary files

| File | Contents |
| --- | --- |
| `index.html` | Ligand, receptor, best score, run ID, Rule-of-Five result, and QC flags |
| `01_input_provenance.html` | Input SMILES, canonical identity, run ID, and receptor source |
| `02_ligand_preparation.html` | Conformer generation, Meeko preparation, and Rule-of-Five data |
| `03_receptor_preparation.html` | Receptor source, chain selection, repair, preparation, and versions |
| `03a_receptor_selection.html` | Ranked targets, preflight candidates, selected PDB, score, and method |
| `04_docking_results.html` | Vina version, box, replica results, selected replica, and ranked poses |
| `05_distance_contacts.html` | Heavy-atom distance contacts for the selected pose |
| `06_methods_qc.html` | QC flags, software versions, seed, replicas, and pH |
| `07_md_handoff_readiness.html` | Handoff validation and topology-readiness limits |
| `08_raw_data.html` | Complete run state rendered as readable JSON |
| `09_visualization.html` | Interactive receptor-backbone and ligand-pose view plus pose table |
| `run_summary.json` | Compact run ID, ligand, receptor, score, box source, QC, and MD status |

The visualization embeds its coordinates and rendering code. It opens directly from disk in a
modern browser and does not fetch scripts or neighboring structure files. Drag to rotate,
scroll to zoom, and use Reset view to restore the starting orientation.

Package creation checks the HTML5 structure, local links, template substitution, summary
fields, and embedded viewer data before publication.

### Overlay files

| File | Contents |
| --- | --- |
| `all_best_poses.sdf` | Chemistry-authoritative selected ligand pose |
| `receptor_overlay_reference.pdb` | Exact prepared receptor used for this run |
| `ligand_overlay.json` | Receptor profile, run identity, score, pose path, and display color |

### `MD_Handoff` files

| Group | Files |
| --- | --- |
| Receptor | `receptor_source.cif`, `prepared_receptor.pdb`, `receptor.pdbqt`, `receptor.json`, `protein_stats.json`, `receptor_prep.log` |
| Ligand input | `ligand_input.sdf`, `ligand.pdbqt`, `ligand_2d.svg` |
| Docked poses | `docked.pdbqt`, `docked_poses.sdf`, `docked_replica_<n>.pdbqt`, `best_pose_ligand.sdf`, `best_pose_ligand.pdb`, `complex_best_pose.pdb` |
| Results | `docking.log`, `vina_replica_<n>.log`, `docking_results.json`, `interactions.json`, `meeko_export.log` |
| Provenance | `atom_mapping.json`, `provenance.json`, `README_MD_HANDOFF.md`, `manifest.json` |

Use `best_pose_ligand.sdf` when chemical bond order matters. `manifest.json` records file sizes
and SHA-256 checksums. The complete contract is in [`../docs/OUTPUTS.md`](../docs/OUTPUTS.md).

## Resume and cleanup

Nextflow keeps `work/` so completed tasks can be reused:

```bash
nextflow run nextflow/main.nf -profile local \
  --input nextflow/examples/input.smi -resume
```

Delete cached work only after it is no longer needed:

```bash
nextflow clean -f
```

`work/`, `results/`, `.nextflow/`, and `.nextflow.log*` are ignored by Git.

## Validation

Run the DAG without scientific software execution:

```bash
nextflow run nextflow/main.nf -profile local,test -stub-run \
  --input nextflow/examples/input.test.smi
```

Run repository tests and lint checks from the root:

```bash
python3 -m unittest discover -s tests -v
ruff check nextflow/bin tests
```

## Scientific boundary

The workflow performs ligand and receptor preparation, rigid-receptor docking, ranked-pose
export, a geometric contact screen, reporting, and structural packaging. It does not perform
induced-fit docking, molecular dynamics, free-energy calculations, toxicity prediction, or
experimental validation. Review [`../docs/SCIENTIFIC_METHOD.md`](../docs/SCIENTIFIC_METHOD.md)
before using a result in a study.
