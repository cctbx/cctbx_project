import sys
output = sys.stdout
from scitbx.matrix import dihedral_angle


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


def build_dihedrals(mainchain): 
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


