#        Copyright 2021  Richardson Lab at Duke University
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

from __future__ import nested_scopes, generators, division, absolute_import
from __future__ import  with_statement, print_function
import os, sys

from iotbx.data_manager import DataManager    #   Load in the DataManager
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


class Failure(Exception):
  message = ""

  def __init__(self, msg):
    self.message=msg


def main(inFile):
  # for testing only
  alt = 'A'
  mgr = loadModel(inFile)

  residues = getResidueDihedrals(mgr, alt,
                                 name=os.path.splitext(inFile)[0])
  for r in residues:
    #print(residueString(r))
    pass


def loadModel(filename):
  dm = DataManager()             #   Initialize the DataManager and call it dm
  dm.set_overwrite(True)         #   tell the DataManager to overwrite files with the same name
  #print("Reading file")
  manager = dm.get_model(filename)
  return manager


def getResidueDihedrals(mgr, altcode='', name='', errorFile = sys.stderr):
  # use altcode='A', 'B', etc. for real alt work
  # if name is provided produce a debug file <name>.suites
  global manager, hierarchy, backbone_hierarchy, chain

  if altcode=='':  altcode = 'A'
  manager = mgr
  hierarchy = manager.get_hierarchy()

  # dboutput = open("hierarchy.show.txt", "w")
  selector = "name P or name O5' or name C5' or name C4' or name C3' or name O3'"
  selection = manager.selection(selector)
  #print (len(selection), "atoms")
  backbone_hierarchy = hierarchy.select(selection)
  # backbone_hierarchy.show(dboutput)
  chains = backbone_hierarchy.chains()
  all_residues = []
  cchains = list(chains)
  nchains = len(cchains)
  # print(cchains)
  if nchains == 0:
    print("Model contains no RNA", file=errorFile)
    return []

  for chain in cchains:
    # print("chain ", chain.id)
    conf = get_matching_conformer(chain, altcode)

    names = conf.get_residue_names_padded(pad=False)
    if all((n == "HOH" for n in names)):
      # print("ignoring water chain")
      continue
    # print(list(names))
    ids = conf.get_residue_ids(pad=False)
    # #print(list(ids))
    backbone_atoms = conf.atoms()

    try:
      residues = build_dihedrals(backbone_atoms, chain.id, conf.altloc)
      all_residues.extend(residues)
    except Failure as e:
      print("chain ", chain.id, e.message, file=errorFile)
      continue

  # if name != "":  # useful for seeing what suites were generated
  #   writeSuitesFile(residues, name)
  if len(all_residues) == 0:
    errorFile.write("read no residues: perhaps wrong alternate code\n")
  return all_residues


def get_matching_conformer(chain, altcode):
  confo = [c for c in chain.conformers() if c.altloc in (altcode, '')]
  if len(confo)>0:
    conf = confo[0]
  else:
    conf = chain.only_conformer()
  return conf


def build_dihedrals(mainchain, chain_id, alt):
  residues = []
  # print (list(mainchain.extract_name()))  # gives a list of atom names

  backbone = [atom for atom in mainchain if atom.name in bone]
  residue = None

  k = 0
  i = find_start(backbone)
  while k < len(backbone) - 3:
    pivot = backbone[k + 1]
    labels = pivot.fetch_labels()
    res_number = int(labels.resseq)
    if residue is None or res_number != residue.sequence:
      # we are seeing our first angle of a new residue
      residue = new_residue(residues, pivot)
    name, angle = easy_make_dihedral(backbone, i, k)
    # ----------------- on failure: new strategy ----------------------
    if angle == 9999:
      # this residue is strangely ordered due to partial alt conformers;
      # build dihedrals the hard way hereafter
      k = harder_build_chain(backbone, i, k, residue, residues,
          chain_id, alt)
      return residues

    else:
      residue.angle[i] = angle
      k = k + 1
      i = (i + 1) % 6
  return residues


def find_start(backbone):
  for i in range(len(bones)):
    if backbone[0].name == bones[i]:
      return i
  assert False, "backbone atom not recognized"


def easy_make_dihedral(backbone, i, k):
  dh = dihedrals[i % 6]
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


def harder_build_chain(backbone, i, k, residue, residues, chain_id, alt):
  if k - i - 2 < 0:
    syn_chain, i = jump_start()
    k = 0
  else:
    syn_chain = backbone[k:k + 3]
  ksyn = k  # point on backbone where syn_chain shadowing begins
  # print("syn_chain")
  # for atom in syn_chain:
  #   print(atom.name, atom.fetch_labels().resseq, atom.serial) #just for debugging

  ii = i  # indicates which dihedral angle we are working on
  while k < len(backbone) - 3:
    k2 = harder_build_residue(syn_chain, ii, k-ksyn, residue)
    # k2 is an index into syn_chain
    if k2 > k-ksyn:
      k = k2+ksyn
      ii = 0
      residue = new_residue(residues, syn_chain[k2+1])
    else:
      # things are horribly broken, just keep on looking for sanity
      ok = resync(residues, syn_chain, k2)
      if not ok:  # we can't fix it, terminate here
        break
      k += 1
  # delete a possible final dead one
  if residues[-1].is_dead() or \
      (len(residues) > 1 and residues[-1].sequence == residues[-2].sequence):
    del residues[-1]
  return k


def harder_build_residue(syn_chain, i, kk, residue):
  # kk is an index to syn_chain
  assert len(syn_chain) >= kk + 3, "We have run out of syn_chain"

  res = syn_chain[-1].fetch_labels().resseq
  res_atoms, res_names = get_one_residue(res)
  res2_atoms, res2_names = get_one_residue(str(int(res)+1))

  for j in range(i, 6):
    atoms = res_atoms if j < 4 else res2_atoms
    names = res_names if j < 4 else res2_names
    if bones[j+3] in names:  # add atom 3 ahead of dihedral start
      # so that we have 4 atoms to work with
      x = names.index(bones[j+3])
      atom = atoms[x]
      syn_chain.append(atom)
      # print("appending", atom.name, atom.fetch_labels().resseq, atom.serial)
      name, angle = easy_make_dihedral(syn_chain, j, kk)
      # print(" ", syn_chain[kk+1].fetch_labels().resseq, name, angle, kk, j)
      residue.angle[j] = angle  # if second failure, we get 9999
      kk = kk + 1
    # syn_name = [a.name for a in syn_chain]
  assert len(syn_chain) >= kk + 3, "We have run out of syn_chain"
  return kk


def resync(residues, syn_chain, k2):
  "previous residue is not working for us, find next useful one"
  oldres = syn_chain[-1].fetch_labels().resseq
  for r in range(1, 21):
    res = str(int(oldres) + r)
    selection=manager.selection("resseq " + res + " and (name P or name O5' or name C5' or name C4' or name C3' or name O3')")
    res_hierarchy = hierarchy.select(selection)
    atoms = res_hierarchy.atoms()
    if len(atoms) > 0:  # yes! we have found a useful residue
      syn_chain.append(atoms[0])
      return True
  return False


def jump_start():
  """Used when a hiccup occurs at the beginning of the chain.
  Get three items on the start of syn_chain so that harder_build_residue
  has an adequate starting point to work with.
  """
  syn_chain = []
  res = list(hierarchy.residue_groups())[0].resseq
  res_atoms, res_names = get_one_residue(res)
  res2_atoms, res2_names = get_one_residue(str(int(res)+1))

  i=1
  while i<7:
    if bones[i] in res_names:
      # we have found our starting point
      atom = res_atoms[res_names.index(bones[i])]
      syn_chain.append(atom)
      break
    i += 1
  j = i+1; n = 1
  jmax = 10  # if we haven't got them by now, give up
  while n<3 and j<jmax:
    atoms = res_atoms if j < 6 else res2_atoms
    names = res_names if j < 6 else res2_names
    if bones[j] in names:  # add atom 3 ahead of residue start
      # so that we have 4 atoms to work with
      atom = atoms[names.index(bones[j])]
      syn_chain.append(atom)
      n += 1
    j += 1
  if j>=jmax:
    raise Failure(" unable to find a backbone")
  return syn_chain, i


def new_residue(residues, pivot):
  """Start a new Residue object, whose angles will be filled in later as we
     continue down the chain."""
  residue = Residue("", "", [9999, 9999, 9999, 9999, 9999, 9999])
  labels = pivot.fetch_labels()
  res_number = int(labels.resseq)
  alt = labels.altloc
  chainID = labels.chain_id
  residue.sequence = res_number
  grandpa = pivot.parent().parent()
  residue_name = grandpa.unique_resnames()[0].rjust(3)
  # CIF files sometimes give shorter names, but PDB files always give 3 chars
  # rjust bridges the gap
  id = (" ", "1", "{:>2}".format(chainID), labels.resseq, " ", "{:1}".format(alt), residue_name)
  residue.pointIDs = id
  base = findBase(residue_name)
  residue.base = base
  if residue.base is not None:
    residues.append(residue)
  return residue


def get_one_residue(res):
  selection=manager.selection("resseq " + res)
  res_hierarchy = hierarchy.select(selection)
  atoms = res_hierarchy.atoms()
  names = [a.name for a in atoms]
  # print("residue ", res)
  # for a in atoms:
  #   print("  ", a.name, a.parent().parent().resseq)
  return atoms, names


# ----------------------------  tool room  -----------------------------

# useful for seeing what suites were generated
def writeSuitesFile(r, name):
  file = open(name + ".suites", "w")
  for r in all_residues:
    id = ":".join(r.pointIDs)
    angles = "".join([":{:.3f}".format(a) for a in r.angle])
    # print(id + angles, file=file)
  file.close()


# for debug printouts
def nameList(atoms):
  list = [a.name for a in atoms]
  return list


def residueString(r):
  sizes=[1, 1, 1, 3, 1, 1, 3]
  id = "".join(["{num:{width}}:".format(num=x,width=y) for x, y in zip(r.pointIDs, sizes)])
  angles = "".join(["{:8.3f}:".format(a) for a in r.angle])
  return id + angles[:-1]


# The following code was used to build the dihedrals list above
def make_groups(output):
  i = 0
  for name in names:
    output.write('("{}" ('.format(name))
    for j in range(4):
      output.write('"{}".format(a[i+j]), ')
    output.write("))\n")
    i = i+1


# The following is for debugging and not for production use:
if __name__ == '__main__':
  main(sys.argv[1])
