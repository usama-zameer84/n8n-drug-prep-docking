#!/usr/bin/env python3
"""Download, repair, and prepare the selected receptor."""

import argparse
import hashlib
import json
import os
import subprocess
import sys
import urllib.parse
import urllib.request


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--params", required=True)
    args = ap.parse_args()
    d = json.load(open(args.params))
    from pdbfixer import PDBFixer
    from openmm.app import PDBFile
    from openmm import unit
    from rdkit import Chem, DataStructs
    import meeko
    import rdkit
    import openmm

    pdb_id = d["pdb_id"]
    source_cif = "receptor_source.cif"
    url = f"https://files.rcsb.org/download/{pdb_id}.cif"
    with urllib.request.urlopen(url, timeout=120) as r:
        source_bytes = r.read()
    if len(source_bytes) < 500:
        raise SystemExit("RCSB returned an unexpectedly small mmCIF file")
    open(source_cif, "wb").write(source_bytes)

    fixer = PDBFixer(filename=source_cif)
    chains = list(fixer.topology.chains())
    available = [c.id or str(i) for i, c in enumerate(chains)]
    protein_names = set(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL ASX GLX SEC PYL DA DC DG DT A C G U".split()
    )
    water_names = set("HOH WAT H2O DOD".split())
    heterogens, waters = [], []
    for residue in fixer.topology.residues():
        item = {"name": residue.name, "chain": residue.chain.id, "id": residue.id}
        if residue.name in water_names:
            waters.append(item)
        elif residue.name not in protein_names:
            heterogens.append(item)

    # binding-site reference: highest RDKit Tanimoto similarity among eligible co-crystallized components
    query_mol = Chem.MolFromSmiles(d["smiles"])
    query_fp = Chem.RDKFingerprint(query_mol) if query_mol is not None else None
    excluded = set(
        "HOH WAT H2O DOD NA CL K CA MG MN ZN FE CU CO NI CD HG SO4 PO4 NO3 NH4 GOL EDO PEG PG4 DMS ACT ACE".split()
    )
    candidates = []
    for residue in fixer.topology.residues():
        if (
            residue.name in protein_names
            or residue.name in water_names
            or residue.name in excluded
        ):
            continue
        heavy = [
            a
            for a in residue.atoms()
            if a.element is not None and a.element.symbol != "H"
        ]
        if len(heavy) < 5:
            continue
        coords = []
        for atom in heavy:
            p = fixer.positions[atom.index].value_in_unit(unit.angstrom)
            coords.append([float(p.x), float(p.y), float(p.z)])
        sim, csmiles, cname = None, None, residue.name
        try:
            cu = f"https://data.rcsb.org/rest/v1/core/chemcomp/{urllib.parse.quote(residue.name)}"
            req = urllib.request.Request(
                cu,
                headers={
                    "Accept": "application/json",
                    "User-Agent": "nextflow-drug-prep-docking/1.0 pocket-selection",
                },
            )
            with urllib.request.urlopen(req, timeout=25) as r:
                cd = json.load(r)
            desc = cd.get("rcsb_chem_comp_descriptor") or {}
            csmiles = desc.get("SMILES_stereo") or desc.get("SMILES")
            cname = (cd.get("chem_comp") or {}).get("name") or residue.name
            cmol = Chem.MolFromSmiles(csmiles) if csmiles else None
            if query_fp is not None and cmol is not None:
                sim = float(
                    DataStructs.TanimotoSimilarity(query_fp, Chem.RDKFingerprint(cmol))
                )
        except Exception:
            pass
        center = [sum(row[i] for row in coords) / len(coords) for i in range(3)]
        candidates.append(
            {
                "component_id": residue.name,
                "component_name": cname,
                "chain": residue.chain.id,
                "residue_id": residue.id,
                "heavy_atom_count": len(heavy),
                "query_similarity": round(sim, 4) if sim is not None else None,
                "center_A": [round(v, 3) for v in center],
                "component_smiles": csmiles,
            }
        )
    candidates.sort(
        key=lambda x: (
            -(x["query_similarity"] if x["query_similarity"] is not None else -1.0),
            -x["heavy_atom_count"],
            x["component_id"],
        )
    )
    binding_site_ref = None
    if candidates:
        binding_site_ref = dict(candidates[0])
        binding_site_ref["selection_method"] = (
            "Highest RDKit fingerprint similarity to eligible co-crystallized non-polymer components."
        )
        binding_site_ref["candidate_count"] = len(candidates)

    if d["chain_ids"]:
        missing = [c for c in d["chain_ids"] if c not in available]
        if missing:
            raise SystemExit(
                "Requested chains not present: "
                + ",".join(missing)
                + "; available: "
                + ",".join(available)
            )
        remove = [
            i for i, c in enumerate(chains) if (c.id or str(i)) not in d["chain_ids"]
        ]
        if remove:
            fixer.removeChains(remove)

    fixer.findMissingResidues()
    missing_res = {str(k): list(v) for k, v in fixer.missingResidues.items()}
    if not d["add_missing_residues"]:
        fixer.missingResidues = {}
    fixer.findNonstandardResidues()
    fixer.replaceNonstandardResidues()
    fixer.removeHeterogens(keepWater=d["heterogen_policy"] == "keep_water")
    fixer.findMissingAtoms()
    missing_atom_count = sum(len(atoms) for atoms in fixer.missingAtoms.values())
    missing_atom_count += sum(len(atoms) for atoms in fixer.missingTerminals.values())
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(float(d["ph"]))

    prepared_pdb = "prepared_receptor.pdb"
    with open(prepared_pdb, "w") as fh:
        PDBFile.writeFile(fixer.topology, fixer.positions, fh, keepIds=True)
    receptor_pdbqt = "receptor.pdbqt"
    receptor_json = "receptor.json"
    cmd = [
        sys.executable,
        "-m",
        "meeko.cli.mk_prepare_receptor",
        "--read_pdb",
        prepared_pdb,
        "--write_pdbqt",
        receptor_pdbqt,
        "--write_json",
        receptor_json,
        "--charge_model",
        "gasteiger",
        "--allow_bad_res",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=300)
    open("receptor_prep.log", "w").write(
        (proc.stdout or "") + "\n" + (proc.stderr or "")
    )
    if (
        proc.returncode != 0
        or not os.path.exists(receptor_pdbqt)
        or os.path.getsize(receptor_pdbqt) == 0
    ):
        raise SystemExit("Meeko receptor preparation failed")

    atom_count, residue_keys = 0, set()
    with open(prepared_pdb) as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                atom_count += 1
                residue_keys.add((line[21:22], line[22:27], line[17:20]))

    stats = {
        "source": {
            "pdb_id": pdb_id,
            "format": "PDBx/mmCIF",
            "url": url,
            "sha256": sha256(source_cif),
        },
        "available_chains": available,
        "selected_chains": d["chain_ids"] or available,
        "heterogen_policy": d["heterogen_policy"],
        "non_water_heterogens_detected": heterogens,
        "waters_detected": len(waters),
        "missing_residues_detected": missing_res,
        "missing_residues_built": d["add_missing_residues"],
        "missing_atoms_added": missing_atom_count,
        "prepared_atom_count": atom_count,
        "prepared_residue_count": len(residue_keys),
        "ph": d["ph"],
        "receptor_pdbqt_method": "Meeko mk_prepare_receptor --read_pdb",
        "charge_model": "gasteiger",
        "binding_site_reference": binding_site_ref,
        "binding_site_candidates": candidates,
        "versions": {
            "meeko": str(getattr(meeko, "__version__", "unknown")),
            "rdkit": str(getattr(rdkit, "__version__", "unknown")),
            "openmm": str(getattr(openmm, "__version__", "unknown")),
        },
    }
    json.dump(stats, open("protein_stats.json", "w"), indent=2)

    qc = [
        f
        for f in d.get("qc_flags", [])
        if f != "GRID_CENTER_TO_BE_SELECTED_OR_INFERRED"
    ]
    qc.append("MEEKO_ALLOW_BAD_RES_USED_REQUIRES_REVIEW")
    if d["cx"] is None:
        if binding_site_ref:
            qc.append("GRID_CENTER_FROM_COCRYSTALLIZED_LIGAND_REQUIRES_REVIEW")
        else:
            qc.append("AUTO_CENTER_NOT_SCIENTIFICALLY_VALIDATED")
    if heterogens:
        qc.append("NON_WATER_HETEROGENS_REMOVED_EXPLICITLY")
    if missing_res and not d["add_missing_residues"]:
        qc.append("MISSING_RESIDUES_NOT_MODELLED")
    d["qc_flags"] = sorted(set(qc))
    d["protein_stats"] = stats
    d["auto_grid_center"] = (
        binding_site_ref["center_A"] if d["cx"] is None and binding_site_ref else None
    )
    d["grid_center_source"] = (
        "user_supplied"
        if d["cx"] is not None
        else (
            "co_crystallized_ligand:" + binding_site_ref["component_id"]
            if binding_site_ref
            else "receptor_centroid_fallback"
        )
    )
    json.dump(d, open(args.params, "w"), indent=2)


if __name__ == "__main__":
    main()
