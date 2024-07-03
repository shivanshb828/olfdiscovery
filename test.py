import rdkit.Chem as Chem
from dockstring import load_target
from utils import pdb_to_smiles
import os

smiles = pdb_to_smiles('propionate.pdb')

target = load_target('DRD')
score, aux = target.dock(smiles)

print(aux)

os.environ['PATH'] += ':/Applications/PyMOL.app/Contents/bin'

target.view([aux['ligand']])