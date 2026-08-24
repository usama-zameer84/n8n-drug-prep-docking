# Ligand Docking Workbench

Molecular docking asks a practical question: how might a small molecule sit inside a protein
pocket? A docking program tries different poses and ranks them by score. That makes docking
useful for comparing possible binding modes and deciding what to examine next, but a good score
does not prove that binding happens in the lab.

Getting to that score takes more work than running one command. The ligand and receptor need to
be checked and prepared, file formats have to match, the search box must be defined, and every
result needs enough context to be reviewed later. Ligand Docking Workbench connects those steps.
It prepares one ligand, selects or accepts a receptor, runs AutoDock Vina, and collects the
poses, checks, reports, and structural MD handoff in one reproducible run.

There are two entry points:

| Runner | Best suited to | Entry point | Setup guide |
| --- | --- | --- | --- |
| n8n | Google Drive triggered runs on a self-hosted n8n instance | `n8n/workflow/drug-prep-docking.workflow.json` | [`n8n/README.md`](n8n/README.md) |
| Nextflow | Command-line runs on a workstation or shared-filesystem cluster | `nextflow/main.nf` | [`nextflow/README.md`](nextflow/README.md) |

Each run handles one ligand. `replicas` repeats docking for that ligand; it does not submit
several ligands.

The n8n workflow requires a separate Python task runner. A self-hosted n8n service by itself
cannot execute the workflow's native Python Code nodes. The runner must have the scientific
environment installed and share `/md_project/data` with the n8n service. See the
[`n8n` setup guide](n8n/README.md#requirements) for the complete runner requirements.

This is a research workflow. Docking scores rank predicted poses. They do not establish
binding, efficacy, safety, or clinical usefulness. The MD handoff contains structures and
provenance, but no topology, force-field parameters, solvent, ions, equilibration, or
trajectory.

## Repository layout

```text
n8n/                       n8n workflow, exported Code-node sources, and input examples
nextflow/                  Nextflow workflow, processes, Python helpers, and HTML templates
docs/OUTPUTS.md            complete report and MD-handoff file reference
docs/SCIENTIFIC_METHOD.md  receptor selection, preparation, docking, and interpretation
tests/                     workflow, helper, report, and publication-hygiene tests
environment.yml            shared conda environment
```

Generated runs are written under `results/`, `work/`, or `/md_project/data` and are not part of
the source repository.

## Install the scientific environment

The shared conda environment contains Python 3.12, RDKit, Meeko, PDBFixer, OpenMM, Gemmi,
AutoDock Vina 1.2.7, Nextflow, Java 17, Ruff, and GNU `sha256sum`.

```bash
conda env create -f environment.yml -n drugprep
conda activate drugprep

python -c "import rdkit, openmm, pdbfixer, meeko, gemmi, vina; print('scientific stack OK')"
vina --version
nextflow -version
```

If a shell alias overrides the conda Python, run `unalias python` and check `which python`.

For n8n, install this environment on the Python task-runner host or in its container. Installing
the environment only on the n8n service is not sufficient.

## Input formats

### Nextflow input

Nextflow reads one non-comment record from a `.smi` or text file:

```text
SMILES optional_ligand_name
```

Example:

```text
Cn1c(=O)c2c(ncn2C)n(C)c1=O caffeine
```

Blank lines and lines beginning with `#` are ignored. A second ligand record causes the run to
stop. Set receptor and docking options with Nextflow command-line parameters:

```bash
nextflow run nextflow/main.nf -profile local \
  --input nextflow/examples/input.smi \
  --receptor_selection_mode provided \
  --pdb_id 9H37 \
  --chain_ids A
```

See [`nextflow/README.md`](nextflow/README.md) for profile selection, custom docking boxes, and
resume behavior.

### n8n input

n8n accepts `.txt`, `.smi`, `.smiles`, and `.json` files from the configured Google Drive input
folder. Text files use the first non-empty, non-comment record in the same format shown above.
Text input uses automatic receptor selection and default settings.

JSON input accepts the full configuration. Only `smiles` is required:

```json
{
  "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O",
  "ligand_name": "caffeine",
  "receptor_selection_mode": "provided",
  "pdb_id": "9H37",
  "chain_ids": ["A"],
  "size_x": 20,
  "size_y": 20,
  "size_z": 20,
  "ph": 7.4,
  "exhaustiveness": 8,
  "num_modes": 9,
  "replicas": 1,
  "seed": 20260826
}
```

Ready-to-use inputs are in [`n8n/examples/`](n8n/examples/).

## Input fields

The two runners use the same scientific settings. n8n reads them from JSON. Nextflow reads the
SMILES and ligand name from its input file and takes the remaining values as `--name value`
parameters.

| Field | Default | Accepted value | Purpose |
| --- | --- | --- | --- |
| `smiles` | required | non-empty SMILES, at most 1000 characters | Ligand structure |
| `ligand_name` | `ligand` | text; normalized for filenames | Report and output folder name |
| `receptor_selection_mode` | `auto` | `auto` or `provided` | Select a receptor or use `pdb_id` |
| `pdb_id` | empty | four-character PDB ID | Required in `provided` mode |
| `chain_ids` | all mapped chains | JSON list or comma-separated IDs | Receptor chains to retain |
| `target_organism` | `Homo sapiens` | non-empty name, at most 100 characters | Organism filter for target evidence |
| `target_similarity_threshold` | `70` | 40 to 100 | Minimum ChEMBL similarity percentage |
| `target_candidate_limit` | `5` | integer from 1 to 10 | Maximum ranked target candidates |
| `heterogen_policy` | `remove_all` | `remove_all` or `keep_water` | Receptor heterogen handling |
| `add_missing_residues` | `false` | boolean | Allow PDBFixer to build missing residues |
| `allow_multicomponent` | `false` | boolean | Permit dot-separated components in one SMILES |
| `center_x`, `center_y`, `center_z` | inferred | all three numbers from -10000 to 10000 Å | Explicit docking-box center |
| `size_x`, `size_y`, `size_z` | `20` | 8 to 30 Å each; volume at most 27000 Å³ | Docking-box dimensions |
| `ph` | `7.4` | 4 to 10 | Receptor hydrogen-addition pH |
| `exhaustiveness` | `8` | integer from 8 to 64 | Vina search effort |
| `num_modes` | `9` | integer from 1 to 20 | Maximum poses retained per replica |
| `energy_range` | `3` | 1 to 10 kcal/mol | Vina score window |
| `seed` | `20260824` | integer from 1 to 2147483647 | Base random seed |
| `cpu` | `1` | integer from 1 to 8 | Vina CPU threads |
| `replicas` | `1` | integer from 1 to 3 | Docking repeats for the same ligand |
| `timeout_seconds` | `900` | integer from 60 to 1800 | Timeout for each Vina replica |
| `cutoff` | `4.5` | 2.5 to 6 Å | Heavy-atom contact cutoff |

Supply all three box-center coordinates or omit all three. When no center is supplied, the
workflow uses the most similar eligible co-crystallized component. If none is suitable, it
uses the receptor centroid and records a QC flag.

A dot-separated SMILES still represents one input record. Review salts, protonation,
tautomerism, stereochemistry, metals, cofactors, retained waters, and covalent chemistry before
interpreting a run.

## Processing stages

1. Validate the ligand and run parameters.
2. Canonicalize the ligand, generate a seeded 3D conformer, and calculate descriptors.
3. Select a receptor from ChEMBL and RCSB evidence, or validate the supplied PDB entry.
4. Repair the selected chains with PDBFixer and prepare ligand and receptor PDBQT files with
   Meeko.
5. Select or infer the docking box and run the configured Vina replicas.
6. Parse ranked poses and replica statistics.
7. Calculate a transparent heavy-atom distance-contact screen for the selected pose.
8. Build reports, provenance records, checksums, and the structural MD handoff.

The receptor-selection report records candidate targets, structures, chain mappings,
preparation preflight results, and the final rationale. Details are in
[`docs/SCIENTIFIC_METHOD.md`](docs/SCIENTIFIC_METHOD.md).

## Outputs

Nextflow publishes one run under:

```text
results/<ligand>/
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

n8n uploads the same report package and MD-handoff files to the configured Google Drive report
folder. It also records Drive delivery identifiers in its final summary.

The report pages cover input provenance, ligand preparation, receptor selection and
preparation, docking settings and poses, geometric contacts, quality-control flags, MD
readiness, raw run data, and an interactive structure view. The Nextflow visualization embeds
its coordinates and renderer, so it opens from disk without a server or internet connection.

`MD_Handoff/` contains the source and prepared receptor, ligand input and PDBQT files, docked
poses, the selected ligand in SDF and PDB form, the receptor-ligand complex, docking logs,
interaction data, atom mapping, provenance, and a SHA-256 manifest. Use
`best_pose_ligand.sdf` as the chemistry-authoritative ligand structure.

The complete file-by-file contract is in [`docs/OUTPUTS.md`](docs/OUTPUTS.md).

## Run the checks

```bash
python3 -m unittest discover -s tests -v
ruff check nextflow/bin tests
python3 n8n/tools/extract_code_nodes.py --check
```

The test suite checks workflow structure, exported Code-node consistency, placeholder Drive
configuration, Python and JavaScript syntax, Nextflow process wiring, report completeness,
offline visualization data, and publication hygiene. CI also runs a Nextflow stub execution.

## License and citation

The repository is licensed under MIT. See [`CITATION.cff`](CITATION.cff) for repository citation
metadata. AutoDock Vina, Meeko, RDKit, OpenMM, PDBFixer, Nextflow, n8n, ChEMBL, RCSB, UniProt,
and their data remain subject to their own licenses and citation requirements.
