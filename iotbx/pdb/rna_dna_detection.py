from scitbx import matrix

# geostd values 2009-10-17
bond_distance_ideal_by_bond_atom_names_v3 = {
  "OP1 P":   1.485,
  "OP2 P":   1.485,
  "O5' P":   1.593,
  "C1' C2'": 1.521,
  "C1' O4'": 1.416,
  "C2' O2'": 1.407, # RNA only
  "C2' C3'": 1.529,
  "C3' O3'": 1.422,
  "C3' C4'": 1.528,
  "C4' O4'": 1.456,
  "C4' C5'": 1.511,
  "C5' O5'": 1.427}
bond_distance_ideal_by_bond_atom_names_v2 = {}
for key,value in bond_distance_ideal_by_bond_atom_names_v3.items():
  key = (key
    .replace("OP1", "O1P")
    .replace("OP2", "O2P")
    .replace("'", "*"))
  bond_distance_ideal_by_bond_atom_names_v2[key] = value

# geostd value 2009-10-17 (identical for A, C, G, U)
c1_n_distance_ideal = 1.463

class _make_tables(object):

  def __init__(O):
    O.cache = None

  def __call__(O, index):
    if (O.cache is None):
      O.cache = []
      atoms_t = ("C5'", "C4'", "O4'", "C1'", "C2'", "C3'", "O3'")
      atoms_s = tuple([a.replace("'","*") for a in atoms_t])
      for bonds_1 in [("OP1 P", "OP2 P"), ("O1P P", "O2P P")]:
        bonds_t = set(bonds_1 + (
          "O5' P",
          "C1' C2'", "C2' C3'", "C3' C4'", "C3' O3'", "C4' C5'", "C4' O4'",
          "C1' O4'", "C5' O5'"))
        O.cache.append((atoms_t, bonds_t))
        O.cache.append((atoms_s, set([b.replace("'","*") for b in bonds_t])))
    return O.cache[index]

_tables = _make_tables()

class atom_name_analysis(object):

  __slots__ = [
    "sub_classification",
    "required_atoms",
    "required_bonds",
    "have_o2",
    "c2_o2",
    "c1",
    "have_all_required_atoms"]

  def __init__(O, atom_dict):
    O.sub_classification = None
    O.required_atoms = None
    O.required_bonds = None
    O.have_o2 = None
    O.c2_o2 = None
    O.c1 = None
    O.have_all_required_atoms = False
    if ("P" not in atom_dict): return
    if   ("OP1" in atom_dict):
      if ("OP2" not in atom_dict): return
      tab_offs = 0
    elif ("O1P" in atom_dict):
      if ("O2P" not in atom_dict): return
      tab_offs = 2
    else:
      return
    if   ("O5'" in atom_dict):
      O.required_atoms, O.required_bonds = _tables(index=tab_offs)
      O.have_o2 = "O2'" in atom_dict
      O.c2_o2 = "C2' O2'"
      O.c1 = "C1'"
    elif ("O5*" in atom_dict):
      O.required_atoms, O.required_bonds = _tables(index=tab_offs+1)
      O.have_o2 = "O2*" in atom_dict
      O.c2_o2 = "C2* O2*"
      O.c1 = "C1*"
    else:
      return
    if (tab_offs == 2 and O.c2_o2 == "C2* O2*"):
      O.sub_classification = ""
    elif (tab_offs == 0 and O.c2_o2 == "C2' O2'"):
      O.sub_classification = "v3"
    else:
      return
    for required_atom in O.required_atoms:
      if (not required_atom in atom_dict): return
    O.have_all_required_atoms = True

  def bond_distance_ideal_by_bond_atom_names(O):
    if (O.sub_classification == "v3"):
      return bond_distance_ideal_by_bond_atom_names_v3
    return bond_distance_ideal_by_bond_atom_names_v2

def classification(atom_dict, bond_list):
  O = atom_name_analysis(atom_dict=atom_dict)
  if (not O.have_all_required_atoms):
    return None
  rna_indicator = False
  bonds_matched = set()
  for bond in bond_list:
    pair = [bond.atom_id_1, bond.atom_id_2]
    pair.sort()
    pair = " ".join(pair)
    if (pair in O.required_bonds):
      bonds_matched.add(pair)
    elif (O.have_o2 and pair == O.c2_o2):
      rna_indicator = True
  if (len(bonds_matched) != len(O.required_bonds)):
    return None
  if (rna_indicator): return "RNA"+O.sub_classification
  return "DNA"+O.sub_classification

class residue_analysis(object):

  __slots__ = [
    "atom_dict",
    "atom_name_analysis",
    "long_distances",
    "c1_n_closest",
    "c1_n_closest_distance"]

  def __init__(O, residue_atoms, distance_tolerance=0.5):
    def build_atom_dict():
      result = {}
      for atom in residue_atoms:
        key = atom.name.strip()
        if (key in result):
          return None
        result[key] = atom
      return result
    O.atom_dict = build_atom_dict()
    O.atom_name_analysis = None
    O.long_distances = None
    O.c1_n_closest = None
    O.c1_n_closest_distance = None
    if (O.atom_dict is None):
      return
    O.atom_name_analysis = ana = atom_name_analysis(atom_dict=O.atom_dict)
    if (not ana.have_all_required_atoms):
      return
    O.long_distances = []
    bdibban = ana.bond_distance_ideal_by_bond_atom_names()
    def check_distance(name_pair):
      bond_atom_names = name_pair.split()
      sites = [matrix.col(O.atom_dict[key].xyz) for key in bond_atom_names]
      distance_model = abs(sites[0] - sites[1])
      distance_ideal = bdibban[name_pair]
      if (distance_model > distance_ideal + distance_tolerance):
        O.long_distances.append((name_pair, distance_model))
    for name_pair in ana.required_bonds:
      check_distance(name_pair=name_pair)
    if (ana.have_o2):
      check_distance(name_pair=ana.c2_o2)
    site_c1 = matrix.col(O.atom_dict[ana.c1].xyz)
    O.c1_n_closest_distance = c1_n_distance_ideal + distance_tolerance
    for atom in residue_atoms:
      sw = atom.name.startswith
      if (not (sw("N") or sw(" N"))): continue
      distance = abs(site_c1 - matrix.col(atom.xyz))
      if (    O.c1_n_closest_distance >= distance
          and (   O.c1_n_closest_distance != distance
               or O.c1_n_closest is None)):
        O.c1_n_closest = atom
        O.c1_n_closest_distance = distance

  def is_rna_dna(O):
    if (O.long_distances is None): return False
    if (len(O.long_distances) != 0): return False
    if (O.c1_n_closest is None): return False
    return True

  def is_rna(O):
    if (not O.is_rna_dna()): return False
    if (not O.atom_name_analysis.have_o2): return False
    return True
