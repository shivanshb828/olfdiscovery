import rdkit.Chem as Chem
from dockstring import load_target
from utilities import pdb_to_smiles
import os
from dotenv import load_dotenv

load_dotenv()

smiles = pdb_to_smiles('propionate.pdb')
target = load_target('or51e2', os.getenv('TEMP_DIR'), os.getenv('TARGET_DIR'))
score, aux = target.dock(smiles)
print(aux)

os.environ['PATH'] += ':/Applications/PyMOL.app/Contents/bin'

target.view([aux['ligand']])

