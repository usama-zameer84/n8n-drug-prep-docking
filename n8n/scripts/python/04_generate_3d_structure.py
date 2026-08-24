"""n8n Code node: Generate 3D Structure.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import os
    import subprocess
    import sys
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit.Chem.Draw import rdMolDraw2D

    d = _items[0]["json"]
    out_dir = d["output_dir"]
    mol = Chem.MolFromSmiles(d["smiles"])
    if mol is None:
        raise ValueError("RDKit rejected the SMILES")
    if len(Chem.GetMolFrags(mol)) != 1:
        raise ValueError("Exactly one ligand molecule is required")

    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    inchi_proc = subprocess.run([
        sys.executable, "-c",
        "from rdkit import Chem; import sys; print(Chem.MolToInchiKey(Chem.MolFromSmiles(sys.argv[1])))",
        d["smiles"],
    ], capture_output=True, text=True, timeout=30)
    if inchi_proc.returncode != 0 or not inchi_proc.stdout.strip():
        raise ValueError("RDKit InChIKey generation failed: " + inchi_proc.stderr[-500:])
    inchi_key = inchi_proc.stdout.strip()
    mol_h = Chem.AddHs(mol)
    embed_status = AllChem.EmbedMolecule(
        mol_h,
        randomSeed=int(d["seed"]),
        enforceChirality=True,
        useExpTorsionAnglePrefs=True,
        useBasicKnowledge=True,
        useSmallRingTorsions=True,
        useMacrocycleTorsions=True,
        ETversion=2,
    )
    if embed_status != 0:
        raise ValueError("RDKit ETKDGv3 3D embedding failed")

    if AllChem.MMFFHasAllMoleculeParams(mol_h):
        force_field = "MMFF94s"
        optimize_status = AllChem.MMFFOptimizeMolecule(mol_h, mmffVariant="MMFF94s", maxIters=1000)
    elif AllChem.UFFHasAllMoleculeParams(mol_h):
        force_field = "UFF_FALLBACK"
        optimize_status = AllChem.UFFOptimizeMolecule(mol_h, maxIters=1000)
    else:
        raise ValueError("No MMFF94s or UFF parameters are available for this ligand")
    if optimize_status < 0:
        raise ValueError(force_field + " optimization failed")

    mol_h.SetProp("_Name", d["ligand_name"])
    mol_h.SetProp("canonical_isomeric_smiles", canonical)
    mol_h.SetProp("inchi_key", inchi_key)
    Chem.MolToMolFile(mol_h, os.path.join(out_dir, "ligand_input.sdf"))

    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DSVG(560, 420)
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol)
    drawer.FinishDrawing()
    with open(os.path.join(out_dir, "ligand_2d.svg"), "w") as fh:
        fh.write(drawer.GetDrawingText())

    out = dict(d)
    out["canonical_smiles"] = canonical
    out["inchi_key"] = inchi_key
    out["num_atoms"] = mol_h.GetNumAtoms()
    out["num_heavy_atoms"] = mol_h.GetNumHeavyAtoms()
    out["conformer_generation"] = {
        "method": "RDKit ETKDGv3-equivalent keyword configuration",
        "random_seed": d["seed"],
        "force_field": force_field,
        "optimization_status": optimize_status,
        "optimization_converged": optimize_status == 0,
    }
    if optimize_status > 0:
        out["qc_flags"] = list(out["qc_flags"]) + ["LIGAND_OPTIMIZATION_DID_NOT_FULLY_CONVERGE"]
    return [{"json": out}]
