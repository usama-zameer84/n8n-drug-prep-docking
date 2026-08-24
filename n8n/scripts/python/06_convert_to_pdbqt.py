"""n8n Code node: Convert to PDBQT.

Call run(_items) with the item structure supplied by the n8n Python runner.
"""


def run(_items):
    import os
    import meeko
    from rdkit import Chem
    from meeko import MoleculePreparation
    from meeko import PDBQTWriterLegacy

    d = _items[0]["json"]
    out_dir = d["output_dir"]
    sdf_path = os.path.join(out_dir, "ligand_input.sdf")
    mol = next(iter(Chem.SDMolSupplier(sdf_path, removeHs=False)), None)
    if mol is None:
        raise ValueError("Cannot read the deterministic ligand SDF")
    preparator = MoleculePreparation()
    setups = preparator.prepare(mol)
    if len(setups) != 1:
        raise ValueError("Meeko must produce exactly one ligand setup")
    result = PDBQTWriterLegacy.write_string(setups[0])
    if isinstance(result, tuple):
        pdbqt_text, is_ok, error_message = result
        if not is_ok:
            raise ValueError("Meeko ligand preparation failed: " + str(error_message))
    else:
        pdbqt_text = result
    if "REMARK SMILES" not in pdbqt_text:
        raise ValueError("Meeko PDBQT is missing the SMILES atom-mapping remarks required for lossless pose export")
    pdbqt_path = os.path.join(out_dir, "ligand.pdbqt")
    with open(pdbqt_path, "w") as fh:
        fh.write(pdbqt_text)
    out = dict(d)
    out["ligand_prep"] = {
        "method": "Meeko MoleculePreparation/PDBQTWriterLegacy",
        "meeko_version": str(getattr(meeko, "__version__", "unknown")),
        "charge_model": "Meeko default Gasteiger",
        "smiles_mapping_present": True,
    }
    return [{"json": out}]
