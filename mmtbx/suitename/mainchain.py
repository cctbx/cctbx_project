import sys
output = sys.stdout
from scitbx.matrix import dihedral_angle
from mmtbx.validation import utils
from suitenamedefs import Residue, findBase


n = "alpha beta gamma delta epsilon zeta"
a = "O3' P O5' C5' C4' C3'"

names = n.split(" ")
bones = 2 * a.split(" ")


dihedrals = (
    ("alpha", ("O3'", "P", "O5'", "C5'", )),
    ("beta", ("P", "O5'", "C5'", "C4'", )),
    ("gamma", ("O5'", "C5'", "C4'", "C3'", )),
    ("delta", ("C5'", "C4'", "C3'", "O3'", )),
    ("epsilon", ("C4'", "C3'", "O3'", "P", )),
    ("zeta", ("C3'", "O3'", "P", "O5'", ))
)


def build_dihedrals(mainchain, chainID, alt): 
    residues = []
    # print (list(mainchain.extract_name()))
    # gives a list of atom names

    backbone = [atom for atom in mainchain if atom.name.strip() in bones]
    # sanity = (backbone == backbone2)
    for atom in backbone:
      print(atom.serial, atom.name)  
    i = find_start(backbone)
    residue = None
    for k in range(len(backbone) - 3):
        name, angle = easy_make_dihedral(backbone, i, k)
        pivot = backbone[k + 1]
        labels = pivot.fetch_labels()
        res_number = int(labels.resseq)
        print(res_number, name, angle,)
        if residue is None or res_number != residue.sequence:
            # we are seeing our first of a new residue
            residue = Residue("", "", [9999, 9999, 9999, 9999, 9999, 9999])
            residue.sequence = res_number
            grandpa = pivot.parent().parent()
            # print(grandpa.id_str(),grandpa.resid())
            residue_name = grandpa.unique_resnames()[0]
            print(residue_name)
            id = ("1", chainID, labels.resseq, "", alt, residue_name)
            residue.pointIDs = id
            base = findBase(residue_name)
            residue.base = base
            if base is not None:
                residues.append(residue)
        residue.angle[i] = angle
        i = (i + 1) % 6
    return residues


def find_start(backbone):
    for i in range(len(bones)):
        if backbone[0].name.strip() == bones[i]:
            return i
    assert False, "backbone atom not recognized"
        

def easy_make_dihedral(backbone, i, k):
    dh = dihedrals[i]
    if backbone[k].name.strip() == dh[1][0]:
        group = backbone[k:k+4]
        group_names = tuple((a.name.strip() for a in group))
        if group_names == dh[1]:
            clump = backbone[k:k + 4]
            points = [atom.xyz for atom in clump]
            name = dh[0]
            angle = dihedral_angle(sites=points, deg=True) 
            if angle < 0:  # suitename standard: all angles are positive
                angle += 360.0            
            return name, angle
    return "failure: " + str(group_names), 9999
                
def make_dihedral(backbone, k):
    for dh in dihedrals:
        if backbone[k] == dh[1][0]:
            if len(backbone) >= k + 3 and backbone[k:k + 4] == dh[1]:
                clump = backbone[k:k + 4]
                points = [atom.xyz for atom in clump]
                name = dh[0]
                labels = clump[1].fetch_labels()
                angle = calc_dihedral(points)
                return name, angle
                

# The following code was used to build the dihedrals list above
def make_groups():
    i = 0
    for name in names:
        output.write(f'("{name}" (')
        for j in range(4):
            output.write(f'"{atoms[i+j]}", ')
        output.write("))\n")
        i = i+1


def build_dihedrals1(mainchain): 
    print (mainchain.extract_name())
    backbone = []
    for atom in mainchain:
        if atom.name.strip() in bones:
            backbone.append(atom)
    #help(mainchain)
    backbone2 = [atom for atom in mainchain if atom.name.strip() in bones]
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

