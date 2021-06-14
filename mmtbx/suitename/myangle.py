#                Copyright 2021  Richardson Lab
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import os, sys

from iotbx.data_manager import DataManager    #   Load in the DataManager
from mmtbx.validation import utils
from scitbx.matrix import dihedral_angle
from suitenamedefs import Residue, findBase


n = "alpha beta gamma delta epsilon zeta"
a = "O3' P O5' C5' C4' C3'"
bone = (" O3'", " P  ", " O5'", " C5'", " C4'", " C3'")

names = n.split(" ")
bones = 2 * bone

dihedrals = (
    ("alpha", (" O3'", " P  ", " O5'", " C5'", )),
    ("beta", (" P  ", " O5'", " C5'", " C4'", )),
    ("gamma", (" O5'", " C5'", " C4'", " C3'", )),
    ("delta", (" C5'", " C4'", " C3'", " O3'", )),
    ("epsilon", (" C4'", " C3'", " O3'", " P  ", )),
    ("zeta", (" C3'", " O3'", " P  ", " O5'", ))
)


output = sys.stdout  #debug


afile = r"C:\Users\Ken\Desktop\Richardson\molprobity\modules\cctbx_project\mmtbx\suitename\test\4fen.pdb"
# 1q9a for hard alt test

def main(inFile):
  wd = os.getcwd()
  manager = loadModel(inFile)
  residues = getResidueDihedrals(manager)
  for r in residues:
    print(residueString(r))


def loadModel(filename):
  dm = DataManager()             #   Initialize the DataManager and call it dm
  dm.set_overwrite(True)         #   tell the DataManager to overwrite files with the same name
  print("Reading file")
  manager = dm.get_model(filename)
  return manager


def getResidueDihedrals(manager, altcode=''): # use A or B for real alt work
  global backbone_hierarchy, chain
  hierarchy = manager.get_hierarchy()
  i_seq_name_hash = utils.build_name_hash(pdb_hierarchy=hierarchy)

  output = open("hierarchy.show.txt", "w")
  out2 = open("atoms.show.txt", "w")
  model = hierarchy.models()[0]  #!!!
  selector = "name P or name O5' or name C5' or name C4' or name C3' or name O3'"
  selection = manager.selection(selector) 
  print (len(selection))
  backbone_hierarchy = hierarchy.select(selection)
  hierarchy.show(output)
  #backbone_hierarchy.show(output)
  chains = backbone_hierarchy.chains()
  all_residues = []

  
  for chain in chains:
    print("chain ", chain.id)
    conf = get_matching_conformer(chain, altcode)
      
    names = conf.get_residue_names_padded(pad=False)
    if all((n == "HOH" for n in names)):
      print("ignoring water chain")
      continue
    print(list(names))
    ids = conf.get_residue_ids(pad=False)
    print(list(ids))
    backbone_atoms = conf.atoms()

    print(len(backbone_atoms))
    for atom in backbone_atoms:
      print(atom.pdb_label_columns(), file=out2)  
    residues = build_dihedrals(backbone_atoms, chain.id, altcode)
    all_residues.extend(residues)
    print("chain ", chain.id, " complete")
  out2.close()
  print("-- finished --")
  return all_residues


def get_matching_conformer(chain, altcode):
    conformers = list(chain.conformers())
    confo = [c for c in chain.conformers() if c.altloc==altcode]
    if len(confo)>0: 
      conf = confo[0]
    else:
      print("*****no match for ", altcode)
      conf = chain.only_conformer()
    return conf


# for debug printouts
def residueString(r):
    sizes=[1, 1, 3, 1, 1, 3]
    id = "".join([f"{x:{y}}:" for x, y in zip(r.pointIDs, sizes)])
    angles = "".join([f"{a:8.3f}:" for a in r.angle])
    return id + angles[:-1]


def find_start(backbone):
    for i in range(len(bones)):
        if backbone[0].name == bones[i]:
            return i
    assert False, "backbone atom not recognized"
        

def build_dihedrals(mainchain, chainID, alt): 
    residues = []
    # print (list(mainchain.extract_name()))
    # gives a list of atom names

    backbone = [atom for atom in mainchain if atom.name in bone]
    # sanity = (backbone == backbone2)
    for atom in backbone:
      print(atom.serial, atom.name)  
    i = find_start(backbone)
    residue = None
    for k in range(len(backbone) - 3):
        pivot = backbone[k + 1]
        labels = pivot.fetch_labels()
        res_number = int(labels.resseq)
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
        name, angle = easy_make_dihedral(backbone, i, k)
        print(res_number, name, angle)
        # ----------------- can fail! ----------------------
        if angle == 9999:
            harder_build_dihedrals(backbone, i, k, residue, labels)
        residue.angle[i] = angle
        i = (i + 1) % 6
    return residues


def harder_build_dihedrals(backbone, i, k, residue, labels):
    syn_chain = backbone[k - i - 2:k + 2]
    # everything for this residue plus 2 before and 2 after  !!!!?
    res = labels.resseq
    selection=backbone_hierarchy.selection("resseq " + seq)
    res_hierarchy = backbone_hierarchy.select(selection)
    res_atoms = res_hierarchy.atoms()
    res_names = [a.name for a in res_atoms]
    for j in range(i, 6):  #???
      if bones[i + 3] in res_names:
        x = res_names.index(bones[i+3])
        atom = res_atoms[x]
        syn_chain.append(atom)

def easy_make_dihedral(backbone, i, k):
    dh = dihedrals[i]
    if backbone[k].name == dh[1][0]:
        group = backbone[k:k+4]
        group_names = tuple((a.name for a in group))
        if group_names == dh[1]:
            clump = backbone[k:k + 4]
            points = [atom.xyz for atom in clump]
            name = dh[0]
            angle = dihedral_angle(sites=points, deg=True) 
            if angle < 0:  # suitename standard: all angles are positive
                angle += 360.0            
            return name, angle
        return "failure: " + str(group_names), 9999
    return "failure: " + backbone[k].name + "!= " + dh[1][0], 9999
                
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
            output.write(f'"{a[i+j]}", ')
        output.write("))\n")
        i = i+1



if __name__ == '__main__':
  main(sys.argv[1])
