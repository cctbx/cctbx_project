from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function

# This file has been merged into myangle.py

def build_dihedrals1(mainchain):
    print (mainchain.extract_name())
    backbone = []
    for atom in mainchain:
        if atom.name in bones:
            backbone.append(atom)
    #help(mainchain)
    backbone2 = [atom for atom in mainchain if atom.name in bones]
    sanity = (backbone == backbone2)
    for atom in backbone2:
      print(atom.i_seq, atom.name)
    i = find_start(backbone)
    for k in range(len(backbone) - 3):
        name, angle = easy_make_dihedral(backbone, i, k)
        labels = backbone[k + 1].fetch_labels()
        residue = int(labels.resseq)
        print(residue, name, angle,)
        i = (i + 1) % 6

