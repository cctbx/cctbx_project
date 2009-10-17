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

def classification(atom_dict, bond_list):
  if ("P" not in atom_dict): return None
  if   ("OP1" in atom_dict):
    if ("OP2" not in atom_dict): return None
    tab_offs = 0
  elif ("O1P" in atom_dict):
    if ("O2P" not in atom_dict): return None
    tab_offs = 2
  else:
    return None
  if   ("O5'" in atom_dict):
    required_atoms, required_bonds = _tables(index=tab_offs)
    have_o2 = "O2'" in atom_dict
    c2_o2 = "C2' O2'"
  elif ("O5*" in atom_dict):
    required_atoms, required_bonds = _tables(index=tab_offs+1)
    have_o2 = "O2*" in atom_dict
    c2_o2 = "C2* O2*"
  else:
    return None
  if (tab_offs == 2 and c2_o2 == "C2* O2*"):
    sub_classification = ""
  elif (tab_offs == 0 and c2_o2 == "C2' O2'"):
    sub_classification = "v3"
  else:
    return None
  for required_atom in required_atoms:
    if (not required_atom in atom_dict): return None
  rna_indicator = False
  bonds_matched = set()
  for bond in bond_list:
    pair = [bond.atom_id_1, bond.atom_id_2]
    pair.sort()
    pair = " ".join(pair)
    if (pair in required_bonds):
      bonds_matched.add(pair)
    elif (have_o2 and pair == c2_o2):
      rna_indicator = True
  if (len(bonds_matched) != len(required_bonds)):
    return None
  if (rna_indicator): return "RNA"+sub_classification
  return "DNA"+sub_classification
