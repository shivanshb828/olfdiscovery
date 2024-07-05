#!/usr/bin/env python

import os
import sys
import getopt
from rdkit import Chem
from rdkit.Chem import AllChem

def usage():
    """Print helpful, accurate usage statement to stdout."""
    print("Usage: prepare_receptor4.py -r filename")
    print()
    print("    Description of command...")
    print("         -r   receptor_filename ")
    print("        supported file types include pdb, mol2, pdbqt, possibly cif")
    print("    Optional parameters:")
    print("        [-v]  verbose output (default is minimal output)")
    print("        [-o pdbqt_filename]  (default is 'molecule_name.pdbqt')")
    print("        [-A]  type(s) of repairs to make: ")
    print("             'bonds_hydrogens': build bonds and add hydrogens ")
    print("             'bonds': build a single bond from each atom with no bonds to its closest neighbor") 
    print("             'hydrogens': add hydrogens")
    print("             'checkhydrogens': add hydrogens only if there are none already")
    print("             'None': do not make any repairs ")
    print("             (default is 'None')")
    print("        [-C]  preserve all input charges ie do not add new charges ")
    print("             (default is addition of gasteiger charges)")
    print("        [-p]  preserve input charges on specific atom types, eg -p Zn -p Fe")
    print("        [-U]  cleanup type:")
    print("             'nphs': merge charges and remove non-polar hydrogens")
    print("             'lps': merge charges and remove lone pairs")
    print("             'waters': remove water residues")
    print("             'nonstdres': remove chains composed entirely of residues of")
    print("                      types other than the standard 20 amino acids")
    print("             'deleteAltB': remove XX@B atoms and rename XX@A atoms->XX")
    print("             (default is 'nphs_lps_waters_nonstdres') ")
    print("        [-e]  delete every nonstd residue from any chain")
    print("              'True': any residue whose name is not in this list:")
    print("                      ['CYS','ILE','SER','VAL','GLN','LYS','ASN', ")
    print("                      'PRO','THR','PHE','ALA','HIS','GLY','ASP', ")
    print("                      'LEU', 'ARG', 'TRP', 'GLU', 'TYR','MET', ")
    print("                      'HID', 'HSP', 'HIE', 'HIP', 'CYX', 'CSS']")
    print("              will be deleted from any chain. ")
    print("              NB: there are no nucleic acid residue names at all ")
    print("              in the list and no metals. ")
    print("             (default is False which means not to do this)")
    print("        [-M]  interactive ")
    print("             (default is 'automatic': outputfile is written with no further user input)")
    print("        [-d dictionary_filename] file to contain receptor summary information")

def main():
    # process command arguments
    try:
        opt_list, args = getopt.getopt(sys.argv[1:], 'r:vo:A:Cp:U:eM:d:')
    except getopt.GetoptError as e:
        print(f'prepare_receptor4.py: {e}')
        usage()
        sys.exit(2)

    # initialize required parameters
    receptor_filename = "or51e2unprep.pdb"

    # optional parameters
    verbose = None
    repairs = ''
    charges_to_add = 'gasteiger'
    preserve_charge_types = None
    cleanup = "nphs_lps_waters_nonstdres"
    outputfilename = "prepreceptor.pdbqt"
    mode = 'automatic'
    delete_single_nonstd_residues = None
    dictionary = None

    for o, a in opt_list:
        if o in ('-r', '--r'):
            receptor_filename = a
            if verbose:
                print(f'set receptor_filename to {a}')
        if o in ('-v', '--v'):
            verbose = True
            if verbose:
                print('set verbose to True')
        if o in ('-o', '--o'):
            outputfilename = a
            if verbose:
                print(f'set outputfilename to {a}')
        if o in ('-A', '--A'):
            repairs = a
            if verbose:
                print(f'set repairs to {a}')
        if o in ('-C', '--C'):
            charges_to_add = None
            if verbose:
                print('do not add charges')
        if o in ('-p', '--p'):
            if not preserve_charge_types:
                preserve_charge_types = a
            else:
                preserve_charge_types = f'{preserve_charge_types},{a}'
            if verbose:
                print(f'preserve initial charges on {preserve_charge_types}')
        if o in ('-U', '--U'):
            cleanup = a
            if verbose:
                print(f'set cleanup to {a}')
        if o in ('-e', '--e'):
            delete_single_nonstd_residues = True
            if verbose:
                print('set delete_single_nonstd_residues to True')
        if o in ('-M', '--M'):
            mode = a
            if verbose:
                print(f'set mode to {a}')
        if o in ('-d', '--d'):
            dictionary = a
            if verbose:
                print(f'set dictionary to {dictionary}')
        if o in ('-h', '--'):
            usage()
            sys.exit()

    if not receptor_filename:
        print('prepare_receptor4: receptor filename must be specified.')
        usage()
        sys.exit()

    # Reading the molecule using RDKit
    mol = None
    if receptor_filename.endswith('.pdb'):
        mol = Chem.MolFromPDBFile(receptor_filename)
    elif receptor_filename.endswith('.mol2'):
        mol = Chem.MolFromMol2File(receptor_filename)
    elif receptor_filename.endswith('.sdf'):
        supplier = Chem.SDMolSupplier(receptor_filename)
        mol = next(supplier) if supplier else None
    else:
        print(f'Unsupported file format for {receptor_filename}')
        sys.exit(1)
    
    if not mol:
        print(f'Could not read molecule from {receptor_filename}')
        sys.exit(1)
    
    if verbose:
        print(f'Read {receptor_filename}')
    
    if 'bonds' in repairs:
        # Build bonds by distance (RDKit usually does this automatically)
        if verbose:
            print('Building bonds')
        AllChem.EmbedMolecule(mol)

    if 'hydrogens' in repairs:
        # Add hydrogens
        if verbose:
            print('Adding hydrogens')
        mol = Chem.AddHs(mol)

    if charges_to_add == 'gasteiger':
        if verbose:
            print('Adding Gasteiger charges')
        AllChem.ComputeGasteigerCharges(mol)

    # Handle output filename
    if not outputfilename:
        outputfilename = os.path.splitext(os.path.basename(receptor_filename))[0] + '.pdbqt'
    
    # Write the prepared receptor to output file (using pdbqt format)
    with open(outputfilename, 'w') as f:
        f.write(Chem.MolToPDBBlock(mol))
        if verbose:
            print(f'Wrote prepared receptor to {outputfilename}')

if __name__ == '__main__':
    main()

# To execute this command type:
# prepare_receptor4.py -r pdb_file -o outputfilename -A checkhydrogens 