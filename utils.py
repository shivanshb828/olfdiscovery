import rdkit.Chem as Chem

def pdb_to_smiles(pdb_file_path):
    mol = Chem.MolFromPDBFile(pdb_file_path)
    if mol:
        smiles = Chem.MolToSmiles(mol)
        return smiles
    else:
        return None
