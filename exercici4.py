"""
Generate a list of all CA atoms of a given residue type with coordinates.
Parameters: PDB file name, residue type (1-letter or 3-letter code).
""" 

import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList

# --- Argparse Setup ---
parser = argparse.ArgumentParser(description='Generate a list of all CA atoms of a given residue type with coordinates.')
parser.add_argument('pdb_id', help='4-character PDB ID (e.g. 1UBQ)')
parser.add_argument('--restype', default='A', dest='res_type', help='Residue type (1-letter or 3-letter code, default: A)')

args = parser.parse_args()

pdb_id = args.pdb_id
res_type = args.res_type

# 1. Download file automatically (returns the filepath, e.g., './pdb1ubq.ent')
pdbl = PDBList()
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

# 2. Parse the structure directly
parser = PDBParser(QUIET=True)
st = parser.get_structure(pdb_id, file_path)

# --- CODE ----

# We create a dictionary to map single-letter amino acid codes to three-letter codes
ONE_TO_THREE = {
    'A': 'ALA', 'C': 'CYS', 'D': 'ASP', 'E': 'GLU', 'F': 'PHE',
    'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'K': 'LYS', 'L': 'LEU',
    'M': 'MET', 'N': 'ASN', 'P': 'PRO', 'Q': 'GLN', 'R': 'ARG',
    'S': 'SER', 'T': 'THR', 'V': 'VAL', 'W': 'TRP', 'Y': 'TYR'
}


def normalize_res_type(res_input):
    """Convert 1-letter code to 3-letter code if necessary."""
    res_clean = res_input.strip().upper()
    if len(res_clean) == 1:
        return ONE_TO_THREE.get(res_clean, res_clean)
    return res_clean


target_res_type = normalize_res_type(res_type)

selected_ca_atoms = []

# Iterate over all residues in the structure
for residue in st.get_residues():
    # Check if the residue matches the target type
    if residue.get_resname() == target_res_type:
        # Check if the CA atom exists in this residue
        if 'CA' in residue:
            selected_ca_atoms.append(residue['CA'])

# Print results
print(f"C-alpha atoms for residue type '{target_res_type}':\n")
for ca_atom in selected_ca_atoms:
    res = ca_atom.get_parent()
    chain = res.get_parent()
    coords = ca_atom.get_coord()
    
    print(
        f"Chain: {chain.id} | "
        f"Residue: {res.get_resname()}{res.id[1]} | "
        f"Atom: {ca_atom.get_name()} | "
        f"Coords: [{coords[0]:.3f}, {coords[1]:.3f}, {coords[2]:.3f}]"
    )

# example1: python3 exercici4.py 1UBQ --restype A
# example2: python3 exercici4.py 1UBQ --restype ALA