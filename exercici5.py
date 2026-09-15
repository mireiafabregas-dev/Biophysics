"""
Generate a list of backbone peptide bond connections (C-N bonds < cutoff).
Parameters: PDB file name, optional cutoff distance (defaults to 2.5 Å, standard ~1.33 Å).
"""

import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import PDBList

# --- Argparse Setup ---
parser = argparse.ArgumentParser(description='Find backbone C-N contacts within a distance cutoff.')
parser.add_argument('pdb_id', help='4-character PDB ID (e.g. 1UBQ)')
parser.add_argument('--maxdist', type=float, default=2.5, dest='MAXDIST', help='Max distance for a backbone contact (default: 2.5)')

args = parser.parse_args()

pdb_id = args.pdb_id
MAXDIST = args.MAXDIST

# 1. Download file automatically (returns the filepath, e.g., './pdb1ubq.ent')
pdbl = PDBList()
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

# 2. Parse the structure directly
parser = PDBParser(QUIET=True)
st = parser.get_structure(pdb_id, file_path)

# --- CODE ----

cn_atoms = []

for at in st.get_atoms():
    # Strict atom name check for backbone C and N
    if at.id in {'C', 'N'}:
        cn_atoms.append(at)

# Preparing search
nbsearch = NeighborSearch(cn_atoms)

print("NBSEARCH:")

ncontact = 1

for at1, at2 in nbsearch.search_all(MAXDIST):

    if at1 == at2:
        continue

    print(f"Contact: {ncontact}")
    print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
    print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
    print()
    ncontact += 1

# example: python3 exercici5.py 1UBQ --maxdist 2.5

""" 
output example

Contact: 140
at1: <Atom N>, 549, VAL
at2: <Atom C>, 551, VAL

"""
