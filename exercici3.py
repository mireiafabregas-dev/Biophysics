from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList

pdb_id = '1UBQ'  

# 1. Download file automatically (returns the filepath, e.g., './pdb1ubq.ent')
pdbl = PDBList()
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

# 2. Parse the structure directly
parser = PDBParser(QUIET=True)
file_path = pdbl.retrieve_pdb_file('1UBQ', pdir='.', file_format='pdb')
st = parser.get_structure('1UBQ', file_path)

# --- CODE ----

MAXDIST = 3.5  # Define distance for a  contact

# Select only CA atoms
polar_el = {'N','O','S'}

# List to put the polar atoms found
polar_atoms = []

for at in st.get_atoms():

    # at.element contains the chemical symbol ('N', 'O', 'S', 'C', etc.)
    if at.element in polar_el:
        polar_atoms.append(at)
        # print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(polar_atoms)

print("NBSEARCH:")

ncontact = 1

for at1, at2 in nbsearch.search_all(MAXDIST):
    print(f"Contact: {ncontact}")
    print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
    print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
    print()
    ncontact += 1






import argparse
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB import PDBList

# --- Argparse Setup ---
parser = argparse.ArgumentParser(description='Find potential hydrogen bonds and polar contacts.')
parser.add_argument('pdb_id', help='4-character PDB ID (e.g. 1UBQ)')
parser.add_argument('--maxdist', type=float, default=3.5, dest='MAXDIST', help='Max distance for a polar contact (default: 3.5)')

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

# Select only CA atoms
polar_el = {'N','O','S'}

# List to put the polar atoms found
polar_atoms = []

for at in st.get_atoms():

    # at.element contains the chemical symbol ('N', 'O', 'S', 'C', etc.)
    if at.element in polar_el:
        polar_atoms.append(at)
        # print(f"ATOM: {at.get_parent().get_resname()}, {at.get_parent().id[1]}, {at.id}")

# Preparing search
nbsearch = NeighborSearch(polar_atoms)

print("NBSEARCH:")

ncontact = 1

for at1, at2 in nbsearch.search_all(MAXDIST):
    print(f"Contact: {ncontact}")
    print(f"at1: {at1}, {at1.get_serial_number()}, {at1.get_parent().get_resname()}")
    print(f"at2: {at2}, {at2.get_serial_number()}, {at2.get_parent().get_resname()}")
    print()
    ncontact += 1

# example exercici3.py 1UBQ --maxdist 3.5