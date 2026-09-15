"""
Determine the list of pairs of residues whose CA atoms are closer than a given distance
Parameters: PDB file name, distance.
"""

import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList

# --- Argparse Setup ---
parser = argparse.ArgumentParser(description='Find C-alpha contacts under a specified max distance.')
parser.add_argument('pdb_id', help='4-character PDB ID (e.g. 1UBQ)')
parser.add_argument('--maxdist', type=float, default=8.0, dest='MAXDIST', help='Max distance for a contact (default: 8)')

args = parser.parse_args()

pdb_id = args.pdb_id
MAXDIST = args.MAXDIST

# 1. Download file automatically (returns the filepath, e.g., './pdb1ubq.ent')
pdbl = PDBList()
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

# 2. Parse the structure directly
parser = PDBParser(QUIET=True)
st = parser.get_structure(pdb_id, file_path)


# --- CODE ---

select = []

#Select only CA atoms

for at in st.get_atoms():
    if at.id == 'CA':
        select.append(at)
        #print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(select)

print("NBSEARCH:")

#Searching for contacts under HBLNK

ncontact = 1

for at1, at2 in nbsearch.search_all(MAXDIST):
    print(f"Contact: {ncontact}")
    print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
    print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
    print()
    ncontact += 1


# example: python3 exercisi1.py 1UBQ --maxdist 8

""" 
output example

Contact: 326
at1: <Atom CA>, 341, ILE
at2: <Atom CA>, 550, VAL

"""

