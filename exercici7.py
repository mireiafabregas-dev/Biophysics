""" Program to print distances between all atom pairs of two given residues"""

from Bio.PDB.PDBParser import PDBParser
import numpy as np
from Bio.PDB import PDBList, PDBParser

pdb_id = '1UBQ'  

# 1. Download file automatically (returns the filepath, e.g., './pdb1ubq.ent')
pdbl = PDBList()
file_path = pdbl.retrieve_pdb_file(pdb_id, pdir='.', file_format='pdb')

# 2. Parse the structure directly
parser = PDBParser(QUIET=True)
file_path = pdbl.retrieve_pdb_file('1UBQ', pdir='.', file_format='pdb')
st = parser.get_structure('1UBQ', file_path)


# ---- CODE -----


#Selection of Residue 10 of Chain A
res10 = st[0]["A"][10]
#Selection of Residue 20 of Chain A
res20 = st[0]["A"][20]

print("Residue 10 is", res10.get_resname())
print("Residue 20 is", res20.get_resname())

print("\nAtom1 Atom2 dist1 dist2\n-------------------------")
for at10 in res10.get_atoms():      # Replace get_atoms with get_atom if you get an Error!
    for at20 in res20.get_atoms():
        dist = at20 - at10     # Direct procedure with (-) to compute distances
        vector = at20.coord - at10.coord  # Or using numpy coordinates
        distance = np.sqrt(np.sum(vector ** 2))
        print(at10, at20, dist, distance)

center = np.array([10, 10, 10])
print("\nDistance of res10 to {} \n".format(center))
for at10 in res10.get_atoms():
    vect = at10.coord - center
    distance = np.sqrt(np.sum(vect ** 2))
    print(at10, distance)

# example: python3 exercici7.py 1UBQ --res1 10 --res2 20 --chain1 A --chain2 A