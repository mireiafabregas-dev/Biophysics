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

""" 
output example

Residue 10 is GLY
Residue 20 is SER

Atom1 Atom2 dist1 dist2
-------------------------
<Atom N> <Atom N> 26.841194 26.841194
<Atom N> <Atom CA> 27.705828 27.705826
<Atom N> <Atom C> 27.01152 27.011518
<Atom N> <Atom O> 27.704405 27.704405
<Atom N> <Atom CB> 28.923367 28.923365
<Atom N> <Atom OG> 28.661741 28.661741
<Atom CA> <Atom N> 26.76021 26.76021
<Atom CA> <Atom CA> 27.68607 27.68607
<Atom CA> <Atom C> 27.0538 27.0538
<Atom CA> <Atom O> 27.799332 27.799332
<Atom CA> <Atom CB> 28.883638 28.883638
<Atom CA> <Atom OG> 28.57509 28.57509
<Atom C> <Atom N> 26.097872 26.097872
<Atom C> <Atom CA> 27.042477 27.042479
<Atom C> <Atom C> 26.412207 26.412207
<Atom C> <Atom O> 27.175257 27.175257
<Atom C> <Atom CB> 28.193848 28.19385
<Atom C> <Atom OG> 27.822264 27.822264
<Atom O> <Atom N> 26.35525 26.35525
<Atom O> <Atom CA> 27.34624 27.34624
<Atom O> <Atom C> 26.76351 26.76351
<Atom O> <Atom O> 27.56535 27.56535
<Atom O> <Atom CB> 28.469149 28.469149
<Atom O> <Atom OG> 28.049194 28.049194

Distance of res10 to [10 10 10]

<Atom N> 39.46418688865108
<Atom CA> 39.183283863832415
<Atom C> 38.96855544105957
<Atom O> 39.107196877658936


"""
