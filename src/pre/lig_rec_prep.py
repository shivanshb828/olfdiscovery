from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import pybel

def load_and_prepare_molecule(filepath, remove_waters=True, is_receptor=False):
    """
    Load a molecule from a file, optionally remove waters, and prepare by adding hydrogens and optimizing.
    """
    print(f"Loading molecule from: {filepath}")
    if filepath.endswith('.pdb'):
        molecule = Chem.MolFromPDBFile(filepath, removeHs=False)
    elif filepath.endswith('.mol'):
        molecule = Chem.MolFromMolFile(filepath, removeHs=False)
    else:
        raise ValueError("Unsupported file format")

    if is_receptor and remove_waters:
        print("Removing water molecules from receptor...")
        new_molecule = Chem.RWMol(molecule)
        for atom in reversed(new_molecule.GetAtoms()):
            if atom.GetPDBResidueInfo().GetResidueName() == 'HOH':
                new_molecule.RemoveAtom(atom.GetIdx())
        molecule = new_molecule.GetMol()

    print("Adding hydrogens and optimizing molecule...")
    molecule = Chem.AddHs(molecule)
    # print('Computing molecular geometry...')
    # AllChem.EmbedMolecule(molecule)
    # print('Computing force field...')
    # AllChem.MMFFOptimizeMolecule(molecule)
    print("Computing Gasteiger charges...")
    AllChem.ComputeGasteigerCharges(molecule)
    return molecule

def convert_to_pdbqt(molecule, output_filename):
    """
    Convert an RDKit molecule to a PDBQT file using Pybel.
    """
    print(f"Converting molecule to PDBQT and saving to {output_filename}")
    molblock = Chem.MolToMolBlock(molecule)
    pybel_mol = pybel.readstring("mol", molblock)
    # pybel_mol.calccharges(method='gasteiger')
    pybel_mol.write('pdbqt', output_filename, overwrite=True)

def prepare_receptor(filepath):
    """
    Prepare the receptor: load, prepare structure, and save as PDBQT.
    """
    print(f"Preparing receptor from file: {filepath}")
    receptor = load_and_prepare_molecule(filepath, remove_waters=True, is_receptor=True)
    receptor_pdbqt = filepath.replace('.pdb', '_receptor.pdbqt')
    convert_to_pdbqt(receptor, receptor_pdbqt)
    return receptor_pdbqt

def prepare_ligand(filepath):
    """
    Prepare the ligand: load, prepare structure, and save as PDBQT.
    """
    print(f"Preparing ligand from file: {filepath}")
    ligand = load_and_prepare_molecule(filepath, remove_waters=False, is_receptor=False)
    ligand_pdbqt = filepath.replace('.pdb', '_ligand.pdbqt')
    convert_to_pdbqt(ligand, ligand_pdbqt)
    return ligand_pdbqt

# Example usage
receptor_pdbqt = prepare_receptor('or51e2unprep.pdb')
ligand_pdbqt = prepare_ligand('propionate.pdb')

print("Receptor prepared as:", receptor_pdbqt)
print("Ligand prepared as:", ligand_pdbqt)