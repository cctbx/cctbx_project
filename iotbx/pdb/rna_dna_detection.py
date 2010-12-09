from scitbx import matrix
from libtbx import slots_getstate_setstate

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

# 2009-10-17 geostd value for C1'-N, identical for A, C, G, U
c1p_outbound_distance_estimate = 1.463

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

class atom_name_analysis(slots_getstate_setstate):

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

deoxy_ribo_atom_keys = set("C1' C2' C3' O3' C4' O4' C5' O5'".split())

class residue_analysis(slots_getstate_setstate):

  __slots__ = [
    "problems",
    "atom_dict",
    "deoxy_ribo_atom_dict",
    "p_atom",
    "op_atoms",
    "c1p_outbound_candidates",
    "o2p_atom",
    "h_atoms",
    "is_terminus_with_o5p",
    "c1p_outbound_atom",
    "is_rna"]

  def __init__(O, residue_atoms, distance_tolerance=0.5):
    O.problems = []
    O.atom_dict = residue_atoms.build_dict(
      strip_names=True,
      upper_names=True,
      convert_stars_to_primes=True,
      throw_runtime_error_if_duplicate_keys=False)
    if (len(O.atom_dict) != len(residue_atoms)):
      O.problems.append("key_clash")
    deoxy_ribo_problems = []
    O.deoxy_ribo_atom_dict = {}
    for key in deoxy_ribo_atom_keys:
      atom = O.atom_dict.get(key)
      if (atom is None):
        deoxy_ribo_problems.append("missing_"+key)
        continue
      O.deoxy_ribo_atom_dict[key] = atom
      del O.atom_dict[key]
    O.p_atom = None
    for key in O.atom_dict.keys():
      if (key.startswith("P")):
        if (key == "P"):
          O.p_atom = O.atom_dict[key]
          del O.atom_dict[key]
        else:
          O.problems.append("other_P")
    O.op_atoms = []
    for keys in [("OP1", "O1P"),
                 ("OP2", "O2P")]:
      for key in keys:
        atom = O.atom_dict.get(key)
        if (atom is not None):
          O.op_atoms.append(atom)
          del O.atom_dict[key]
          break
      else:
        O.op_atoms.append(None)
        if (O.p_atom is not None):
          O.problems.append("missing_"+keys[0])
    O.o2p_atom = O.atom_dict.get("O2'")
    if (O.o2p_atom is not None):
      del O.atom_dict["O2'"]
    O.h_atoms = {}
    for key in O.atom_dict.keys():
      if (len(key) == 0):
        O.problems.append("blank_name")
        continue
      if (   "HD".find(key[0]) >= 0
          or (    len(key) > 1
              and "HD".find(key[1]) >= 0
              and "0123456789".find(key[0]) >= 0)):
        O.h_atoms[key] = O.atom_dict[key]
        del O.atom_dict[key]
    O.c1p_outbound_candidates = {}
    for key in O.atom_dict.keys():
      if (key.find("'") >= 0):
        if key == "N2'":
          O.deoxy_ribo_atom_dict[key] = atom
          del O.atom_dict[key]
        else:
          O.problems.append("other_prime")
        continue
      if (key.startswith("N") or key.startswith("C")):
        O.c1p_outbound_candidates[key] = O.atom_dict[key]
        del O.atom_dict[key]
    O.is_terminus_with_o5p = False
    if (    len(O.problems) == 0
        and len(deoxy_ribo_problems) != 0
        and len(O.deoxy_ribo_atom_dict) == 1
        and O.deoxy_ribo_atom_dict.get("O5'") is not None):
      O.is_terminus_with_o5p = True
    else:
      O.problems.extend(deoxy_ribo_problems)
    def check_distance(key_pair, site_1, site_2):
      distance_model = abs(site_1 - site_2)
      distance_ideal = bond_distance_ideal_by_bond_atom_names_v3[key_pair]
      if (distance_model > distance_ideal + distance_tolerance):
        return False
      return True
    if (O.p_atom is not None):
      site_1 = matrix.col(O.p_atom.xyz)
      for i,key_pair in enumerate(["OP1 P", "OP2 P"]):
        atom = O.op_atoms[i]
        if (atom is None): continue
        site_2 = matrix.col(atom.xyz)
        if (not check_distance(key_pair, site_1, site_2)):
          O.problems.append("long_distance_P_OP%d" % (i+1))
      atom = O.deoxy_ribo_atom_dict.get("O5'")
      if (atom is not None):
        site_2 = matrix.col(atom.xyz)
        if (not check_distance("O5' P", site_1, site_2)):
          O.problems.append("long_distance_P_O5'")
    if (O.o2p_atom is not None):
      atom = O.deoxy_ribo_atom_dict.get("C2'")
      if (atom is not None):
        site_1 = matrix.col(atom.xyz)
        site_2 = matrix.col(O.o2p_atom.xyz)
        if (not check_distance("C2' O2'", site_1, site_2)):
          O.problems.append("long_distance_C2'_O2'")
    for key_pair in [
          "C1' C2'",
          "C1' O4'",
          "C2' C3'",
          "C3' O3'",
          "C3' C4'",
          "C4' O4'",
          "C4' C5'",
          "C5' O5'"]:
      sites = []
      for key in [key_pair[:3], key_pair[4:]]:
        atom = O.deoxy_ribo_atom_dict.get(key)
        if (atom is None): break
        sites.append(matrix.col(atom.xyz))
      else:
        if (not check_distance(key_pair, *sites)):
          O.problems.append("long_distance_"+key_pair.replace(" ","_"))
    O.c1p_outbound_atom = None
    if (len(O.c1p_outbound_candidates) != 0):
      atom = O.deoxy_ribo_atom_dict.get("C1'")
      if (atom is not None):
        closest_distance = c1p_outbound_distance_estimate + distance_tolerance
        site_1 = matrix.col(atom.xyz)
        for key,atom in O.c1p_outbound_candidates.items():
          site_2 = matrix.col(atom.xyz)
          distance = abs(site_1 - site_2)
          if (    closest_distance >= distance
              and (   closest_distance != distance
                   or O.c1p_outbound_atom is None)):
            O.c1p_outbound_atom = atom
            closest_distance = distance
    if (len(O.problems) == 0):
      O.problems = None
    O.is_rna = (O.o2p_atom is not None)
