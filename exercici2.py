"""
Generate a list of all atoms for a given residue number
Parameters: PDB file name, Residue number (Including Chain if applicable)
"""

import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList

# --- Argparse Setup ---
parser = argparse.ArgumentParser(description='Get all atoms for a given residue number and chain.')
parser.add_argument('pdb_id', help='4-character PDB ID (e.g. 1UBQ)')
parser.add_argument('--resnum', type=int, default=42, dest='RES_NUM', help='Target residue number (default: 42)')
parser.add_argument('--chain', default='A', dest='CHAIN_ID', help='Target chain ID (default: A)')

args = parser.parse_args()

pdb_id = args.pdb_id
RES_NUM = args.RES_NUM
CHAIN_ID = args.CHAIN_ID

# 1. Download file automatically (returns the filepath, e.g., './pdb1ubq.ent')
pdbl = PDBList()
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

# 2. Parse the structure directly
parser = PDBParser(QUIET=True)
st = parser.get_structure(pdb_id, file_path)

# --- CODE ----

selected = []

for at in st.get_atoms():
    res = at.get_parent()
    chain = res.get_parent()
    
    # Match both the integer residue number and chain ID
    if res.id[1] == RES_NUM and chain.id == CHAIN_ID:
        selected.append(at)

print(f"Atoms for Residue {RES_NUM} (Chain {CHAIN_ID}):")
for atom in selected:
    res = atom.get_parent()
    print(f"Residue: {res.get_resname()}{res.id[1]} | Atom: {atom.get_name()} | Coordinates: {atom.get_coord()}")


# example: python3 exercici2.py 1UBQ --resnum 42 --chain A