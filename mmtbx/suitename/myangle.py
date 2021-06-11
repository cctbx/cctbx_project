import sys, os, inspect
import mainchain

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

from iotbx.data_manager import DataManager    #   Load in the DataManager
from mmtbx.validation import utils


afile = r"C:\Users\Ken\Desktop\Richardson\molprobity\modules\cctbx_project\mmtbx\suitename\test\4fen.pdb"
# 1q9a for hard alt test

def main(inFile):
  manager = loadModel(afile)
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
  hierarchy = manager.get_hierarchy()
  i_seq_name_hash = utils.build_name_hash(pdb_hierarchy=hierarchy)

  output = open("hierarchy.show.txt", "w")
  out2 = open("atoms.show.txt", "w")
  model = hierarchy.models()[0]  #!!!
  selector = "name P or name O5' or name C5' or name C4' or name C3' or name O3'"
  selection = manager.selection(selector) 
  print (len(selection))
  backbone_hierarchy = hierarchy.select(selection)
  backbone_hierarchy.show(output)
  chains = backbone_hierarchy.chains()
  all_residues = []
  
  for chain in chains:
    print("chain ", chain.id)
    conformers = list(chain.conformers())
    confo = [c for c in chain.conformers() if c.altloc==altcode]
    if len(confo)>0: 
      conf = confo[0]
    else:
      print("*****no match for ", altcode)
      conf = conformers[0]
      
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
    residues = mainchain.build_dihedrals(backbone_atoms, chain.id, altcode)
    all_residues.extend(residues)
    print("chain ", chain.id, " complete")
  out2.close()
  print("-- finished --")
  return all_residues

# for debug printouts
def residueString(r):
    sizes=[1, 1, 3, 1, 1, 3]
    id = "".join([f"{x:{y}}:" for x, y in zip(r.pointIDs, sizes)])
    angles = "".join([f"{a:8.3f}:" for a in r.angle])
    return id + angles[:-1]


def oldgetResidueDihedrals(model):
  # bb_dihedrals = defaultdict(dict)
  residues = []
  alt_tracker = {}

  pdb_hierarchy = model.get_hierarchy()
  #print(help(pdb_hierarchy))
  grm = model.get_restraints_manager() 
  # print(help(grm))
  # file = open(r"C:\Users\Ken\Desktop\Richardson\temp\model.pdb.geo", "w")
  # header = "# Geometry restraints\n",
  # grm.write_geo_file(
  #       hierarchy = model.get_hierarchy(),
  #       sites_cart=model.get_sites_cart(),
  #       site_labels=model.get_site_labels(),
  #       header=header,
  #       # Stuff for outputting ncs_groups
  #       # excessive_distance_limit = excessive_distance_limit,
  #       # xray_structure=model.get_xray_structure(),
  #       file_descriptor=file)
  # sites_cart = model.get_sites_cart() 
  # #sites_cart = pdb_hierarchy.atoms().extract_xyz()

  # 4/6/2021 office hours ideas:
  selection = model.selection("name *'")  # was "rna backbone"
  print (len(selection))
  # print (selection)
  print (list(selection.as_1d()))  # flex of bool;
  all_atoms = pdb_hierarchy.atoms()
  backbone_atoms = all_atoms.select(selection)
  print(len(backbone_atoms))
  backbone_hierarchy = pdb_hierarchy.select(selection)
  backbone_hierarchy.show()
  for atom in backbone_atoms:
    print(atom.pdb_label_columns)
  
  return

  grm2 = grm.select(selection)
  file2 = open(r"C:\Users\Ken\Desktop\Richardson\temp\model2.pdb.geo", "w")
  grm2.write_geo_file(
        hierarchy = model.get_hierarchy().select(selection),
        sites_cart=model.get_sites_cart().select(selection),
        site_labels=model.get_site_labels().select(selection),
        header="Selected restraints\n",
        # Stuff for outputting ncs_groups
        # excessive_distance_limit = excessive_distance_limit,
        # xray_structure=model.get_xray_structure(),
        file_descriptor=file2)
  dihedral_proxies = grm2.geometry.dihedral_proxies
  print (len(dihedral_proxies))

  # here begins direct copy of utils code:
  i_seq_name_hash = utils.build_name_hash(pdb_hierarchy=pdb_hierarchy)

  def is_blank_or_alt_a(proxy):
    for i in proxy.i_seqs:
       alt = i_seq_name_hash[i][4:5]
       if alt not in [' ', 'A']:
         return False
    return True

  for dp in dihedral_proxies:
    atoms = []
    debug_key = ""
    invert_sign = False
    dp.sort_i_seqs()
    for i in dp.i_seqs:
      atoms.append(i_seq_name_hash[i][0:4].strip())
      debug_key = debug_key+i_seq_name_hash[i]
    if len(atoms) != 4:
      continue
    print(atoms)
    name = utils.match_dihedral_to_name(atoms=atoms)
    #handle dihedral equivalences
    if name == None:
      inverted_atoms = utils.get_inverted_atoms(atoms=atoms, improper=False)
      name = utils.match_dihedral_to_name(atoms=inverted_atoms)
      if name == None:
        inverted_atoms = utils.get_inverted_atoms(atoms=atoms, improper=True)
        name = utils.match_dihedral_to_name(atoms=inverted_atoms)
        if name is not None:
          invert_sign = True
    if (name is not None) and (is_blank_or_alt_a(dp)):
      restraint = geometry_restraints.dihedral(
                                               sites_cart=sites_cart,
                                               proxy=dp)
      key = i_seq_name_hash[dp.i_seqs[1]][4:]
      if alt_tracker.get(key[1:]) is None:
        alt_tracker[key[1:]] = []
      if key[0:1] not in alt_tracker[key[1:]]:
        alt_tracker[key[1:]].append(key[0:1])
      bb_dihedrals[key][name] = restraint.angle_model
      if invert_sign:
        bb_dihedrals[key][name] = bb_dihedrals[key][name] * -1.0
  for key in list(bb_dihedrals.keys()):
    altloc = key[0:1]
    resname = key[1:4]
    chainID = key[4:6]
    resnum = key[6:10]
    i_code = key[10:]
    if 'A' in alt_tracker[key[1:]]:
      if altloc != 'A':
        continue
    if bb_dihedrals[key].get('alpha') is not None:
      alpha = bb_dihedrals[key]['alpha']
    # FIXME will the lookup below ever actually work?
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('alpha') is not None:
      alpha = bb_dihedrals[' '+key[1:]]['alpha']
    else:
      alpha = 9999.
    if bb_dihedrals[key].get('beta') is not None:
      beta = bb_dihedrals[key]['beta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('beta') is not None:
      beta = bb_dihedrals[' '+key[1:]]['beta']
    else:
      beta = 9999.
    if bb_dihedrals[key].get('gamma') is not None:
      gamma = bb_dihedrals[key]['gamma']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('gamma') is not None:
      gamma = bb_dihedrals[' '+key[1:]]['gamma']
    else:
      gamma = 9999.
    if bb_dihedrals[key].get('delta'):
      delta = bb_dihedrals[key]['delta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('delta') is not None:
      delta = bb_dihedrals[' '+key[1:]]['delta']
    else:
      delta = 9999.
    if bb_dihedrals[key].get('epsilon'):
      epsilon = bb_dihedrals[key]['epsilon']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('epsilon') is not None:
      epsilon = bb_dihedrals[' '+key[1:]]['epsilon']
    else:
      epsilon = 9999.
    if bb_dihedrals[key].get('zeta'):
      zeta = bb_dihedrals[key]['zeta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('zeta') is not None:
      zeta = bb_dihedrals[' '+key[1:]]['zeta']
    else:
      zeta = 9999.
    id = ("1", chainID, resnum, i_code, altloc, resname)
    angles = [alpha, beta, gamma, delta, epsilon, zeta]
    for i in range(len(angles)):
      if angles[i] < 0:
        angles[i] += 360.0
    residue = Residue(id, resname, angles)
    residues.append(residue)
  return residues



















def get_rna_backbone_dihedrals(model):
  bb_dihedrals = defaultdict(dict)
  formatted_out = []
  alt_tracker = {}

  pdb_hierarchy = model.get_hierarchy()
  #print(help(pdb_hierarchy))
  geometry = model.get_restraints_manager()  # came up None
  # print(help(geometry))
  sites_cart = model.get_sites_cart() 
  #sites_cart = pdb_hierarchy.atoms().extract_xyz()
  dihedral_proxies = geometry.geometry.dihedral_proxies  #?!!!

  # here begins direct copy of utils code:
  i_seq_name_hash = utils.build_name_hash(pdb_hierarchy=pdb_hierarchy)

  def is_blank_or_alt_a(proxy):
    for i in proxy.i_seqs:
       alt = i_seq_name_hash[i][4:5]
       if alt not in [' ', 'A']:
         return False
    return True

  for dp in dihedral_proxies:
    atoms = []
    debug_key = ""
    invert_sign = False
    dp.sort_i_seqs()
    for i in dp.i_seqs:
      atoms.append(i_seq_name_hash[i][0:4].strip())
      debug_key = debug_key+i_seq_name_hash[i]
    if len(atoms) != 4:
      continue
    name = utils.match_dihedral_to_name(atoms=atoms)
    #handle dihedral equivalences
    if name == None:
      inverted_atoms = utils.get_inverted_atoms(atoms=atoms, improper=False)
      name = utils.match_dihedral_to_name(atoms=inverted_atoms)
      if name == None:
        inverted_atoms = utils.get_inverted_atoms(atoms=atoms, improper=True)
        name = utils.match_dihedral_to_name(atoms=inverted_atoms)
        if name is not None:
          invert_sign = True
    if (name is not None) and (is_blank_or_alt_a(dp)):
      restraint = geometry_restraints.dihedral(
                                               sites_cart=sites_cart,
                                               proxy=dp)
      key = i_seq_name_hash[dp.i_seqs[1]][4:]
      if alt_tracker.get(key[1:]) is None:
        alt_tracker[key[1:]] = []
      if key[0:1] not in alt_tracker[key[1:]]:
        alt_tracker[key[1:]].append(key[0:1])
      bb_dihedrals[key][name] = restraint.angle_model
      if invert_sign:
        bb_dihedrals[key][name] = bb_dihedrals[key][name] * -1.0
  for key in list(bb_dihedrals.keys()):
    altloc = key[0:1]
    resname = key[1:4]
    chainID = key[4:6]
    resnum = key[6:10]
    i_code = key[10:]
    if 'A' in alt_tracker[key[1:]]:
      if altloc != 'A':
        continue
    if bb_dihedrals[key].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[key]['alpha']
    # FIXME will the lookup below ever actually work?
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('alpha') is not None:
      alpha = "%.3f" % bb_dihedrals[' '+key[1:]]['alpha']
    else:
      alpha = '__?__'
    if bb_dihedrals[key].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[key]['beta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('beta') is not None:
      beta = "%.3f" % bb_dihedrals[' '+key[1:]]['beta']
    else:
      beta = '__?__'
    if bb_dihedrals[key].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[key]['gamma']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('gamma') is not None:
      gamma = "%.3f" % bb_dihedrals[' '+key[1:]]['gamma']
    else:
      gamma = '__?__'
    if bb_dihedrals[key].get('delta'):
      delta = "%.3f" % bb_dihedrals[key]['delta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('delta') is not None:
      delta = "%.3f" % bb_dihedrals[' '+key[1:]]['delta']
    else:
      delta = '__?__'
    if bb_dihedrals[key].get('epsilon'):
      epsilon = "%.3f" % bb_dihedrals[key]['epsilon']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('epsilon') is not None:
      epsilon = "%.3f" % bb_dihedrals[' '+key[1:]]['epsilon']
    else:
      epsilon = '__?__'
    if bb_dihedrals[key].get('zeta'):
      zeta = "%.3f" % bb_dihedrals[key]['zeta']
    elif altloc == 'A' and \
         bb_dihedrals[' '+key[1:]].get('zeta') is not None:
      zeta = "%.3f" % bb_dihedrals[' '+key[1:]]['zeta']
    else:
      zeta = '__?__'
    eval = "%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s:%s" \
           % (" ",
              "1",
              chainID,
              resnum,
              i_code,
              altloc,
              resname,
              alpha,
              beta,
              gamma,
              delta,
              epsilon,
              zeta)
    formatted_out.append(eval)
  formatted_out.sort()
  backbone_dihedrals = ""
  for line in formatted_out:
    backbone_dihedrals += line+'\n'
  return backbone_dihedrals


if __name__ == '__main__':
  main(sys.argv[1:])
