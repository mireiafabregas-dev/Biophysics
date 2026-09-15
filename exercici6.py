"""
Generate a list of potential disulfide bonds (S-S contacts between Cys residues).
Parameters: PDB file name, optional cutoff distance (defaults to 3.0 Å).
"""

import argparse
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB import PDBList

# --- Argparse Setup ---
parser = argparse.ArgumentParser(description='Find potential disulfide bonds between Cys residues.')
parser.add_argument('pdb_id', help='4-character PDB ID (e.g. 1UBQ)')
parser.add_argument('--maxdist', type=float, default=3.0, dest='MAXDIST', help='Max distance for a disulfide contact (default: 3.0)')

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

cys_sg_atoms = []

for at in st.get_atoms():
    res = at.get_parent()
    if res.get_resname() == 'CYS' and at.id == 'SG':
        cys_sg_atoms.append(at)


if len(cys_sg_atoms) >= 2:
    # Preparing search
    nbsearch = NeighborSearch(cys_sg_atoms)

    print("NBSEARCH:")

    ncontact = 1

    for at1, at2 in nbsearch.search_all(MAXDIST):
        # Ensure the atoms do not belong to the exact same residue instance
        if at1.get_parent() == at2.get_parent():
            continue
        
        print(f"Contact: {ncontact}")
        print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
        print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
        print()
        ncontact += 1
else:
    print(f"Fewer than 2 Cys SG atoms found in the structure ({len(cys_sg_atoms)} found). No disulfide search performed.")

# example: python3 exercici6.py 4INS --maxdist 3.0

""" 
output example

Contact: 6
at1: <Atom SG>, 464, CYS
at2: <Atom SG>, 497, CYS

"""
