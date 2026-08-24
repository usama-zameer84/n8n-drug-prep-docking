# Scientific method and limits

## Ligand preparation

RDKit parses and canonicalizes the input SMILES, adds hydrogens, generates a deterministic ETKDG conformer, and minimizes it with MMFF where parameters are available. The workflow records the seed, force field, convergence status, atom counts, canonical isomeric SMILES, and InChIKey.

Descriptors include molecular weight, calculated logP, hydrogen-bond donors and acceptors, rotatable bonds, topological polar surface area, and Rule-of-Five violations. These descriptors are filters, not predictions of clinical behavior.

Meeko assigns the ligand PDBQT representation used by Vina. The SDF remains the chemistry-authoritative structure because it preserves element, bond-order, and coordinate information needed for later parameterization.

## Automatic receptor hypothesis

The query ligand is matched against the ChEMBL similarity endpoint. The threshold defaults to 70 percent and can be set from 40 to 100 percent.

For every accepted similar molecule, the workflow retrieves ChEMBL activity records and keeps records that:

- match the requested organism;
- use binding or functional assay types;
- have no data-validity warning;
- resolve to a `SINGLE PROTEIN` target with a UniProt accession.

The target evidence score is:

```text
0.42 * maximum ligand similarity
+ 0.20 * capped activity-record count
+ 0.13 * capped distinct similar-molecule count
+ 0.10 * scaled maximum pChEMBL
+ 0.15 * exact-query evidence indicator
```

The activity-record count is capped at 20. The distinct similar-molecule count is capped at five. The pChEMBL term is scaled between pChEMBL 3 and 9.

The workflow searches RCSB for experimental structures mapped to each target's UniProt accession and containing at least one non-polymer component. The combined structure score is:

```text
0.72 * target evidence score
+ 0.20 * resolution quality
+ 0.08 * non-polymer component indicator
```

Resolution quality scales between 1.5 and 4.5 A and is clamped to zero through one. RCSB polymer-entity references supply the author-chain mapping.

The best-scoring structures are tested in order. A candidate must survive chain selection, PDBFixer processing, and Meeko receptor preparation. The first passing candidate becomes the docking receptor. Failed candidates and their reasons remain in the report.

This ranking favors available assay data and ligand-containing structures. It can reproduce database bias, select an engineered construct, or miss a real target with sparse evidence. Treat the result as a receptor hypothesis.

## Supplied receptor mode

A JSON input can set `receptor_selection_mode` to `provided`, include a four-character `pdb_id`, and optionally specify `chain_ids`. The workflow does not treat a supplied PDB entry as evidence of binding. Receptor identity, construct, assembly, chain, and biological relevance remain the researcher's responsibility.

## Receptor preparation

The workflow downloads the RCSB mmCIF file and records its SHA-256. PDBFixer handles chain selection, nonstandard residues, missing atoms, and hydrogen addition at the configured pH. Building missing residues is off by default because reconstructed loops can distort a binding site.

The default heterogen policy removes waters and non-water components from the prepared receptor. `keep_water` retains waters but still removes other heterogens. Metals, cofactors, structural waters, and catalytic ions require case-specific review.

Meeko creates the receptor PDBQT and receptor JSON files. The workflow records preparation options, versions, atom counts, missing residues, repaired atoms, detected heterogens, and QC flags.

## Binding-site inference

If the user supplies a complete box center, the workflow uses it. Otherwise, it examines eligible co-crystallized non-polymer components. RDKit fingerprints compare each component with the query ligand. The highest Tanimoto similarity wins, with heavy-atom count as the tie breaker. The component centroid becomes the Vina box center.

If no eligible component can be resolved, the receptor centroid is used and a QC flag marks the run as a smoke test. A receptor centroid often misses a real pocket.

## Docking

AutoDock Vina receives the prepared receptor and ligand PDBQT files, box center and dimensions, exhaustiveness, mode count, energy range, CPU count, and deterministic seed. Up to three replicas can run with recorded seeds.

The report stores every selected-replica pose with its Vina score and RMSD bounds. It also stores replica top-score mean, standard deviation, and range. Similar scores do not guarantee similar geometries, so the poses still need visual inspection.

Vina's score supports ranking within a controlled protocol. It is not an experimental affinity and should not be presented as a binding free energy.

## Contact screen

The interaction stage reports receptor atoms and residues within the configured cutoff of the docked ligand. Categories are based on atom identity and distance. The method does not calculate hydrogen-bond geometry, solvation, entropy, polarization, metal coordination, or interaction energy.

## Pose overlays

The overlay registry accepts poses only when the prepared receptor profile matches exactly by source SHA-256, selected chains, heterogen policy, and pH. All accepted poses already share the same coordinate frame. The HTML viewer applies no structural transform.

Overlay appearance can help compare predicted pockets and orientations. It does not show molecular dynamics, occupancy, residence time, or experimental attachment.

## MD boundary

The workflow exports structural inputs for an MD project but does not generate a force-field topology. A safe topology depends on ligand chemistry, protein force field, water model, ions, cofactors, protonation, termini, disulfides, and parameter compatibility.

The handoff status remains `STRUCTURAL_HANDOFF_ONLY_NOT_TOPOLOGY_READY` until those decisions and validation steps are completed outside this workflow.
