import os, sys

import libtbx.load_env

from mmtbx.chemical_components import cif_parser

loaded_cifs = {}

def find_data_dir():
  for relative_path in ["chemical_components",
                        "ext_ref_files/chemical_components"]:
    result = libtbx.env.find_in_repositories(relative_path=relative_path)
    if (result is not None): return result
  return None

data_dir = find_data_dir()

def get_cif_filename(code):
  if (data_dir is None): return ""
  return os.path.join(
    data_dir, "%s" % code[0].lower(), "data_%s.cif" % code.upper())

def is_code(code):
  filename = get_cif_filename(code)
  if os.path.exists(filename): return True
  return False

def is_residue_specified(code, alternate=False):
  if is_code(code):
    return get_atom_names(code, alternate=alternate)
  return False

def get_cif_dictionary(code):
  if code in loaded_cifs:
    cif = loaded_cifs[code]
  else:
    filename = get_cif_filename(code)
    cif = cif_parser.run(filename)
    loaded_cifs[code] = cif
  return cif

def get_alternative_name(code, name):
  for atom, alt in zip(get_atom_names(code),
                       get_atom_names(code, alternate=True)):
    if atom==name: return alt

def get_smiles(code):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  desc = cif.get("_pdbx_chem_comp_descriptor", {})
  for item in desc:
    try:
      if item.program=="CATVS" and item.type=="SMILES_CANONICAL":
        return item.descriptor
    except: pass
  for item in desc:
    try:
      if item.program=="CATVS" and item.type=="SMILES":
        return item.descriptor
    except: pass
  try: return item.descriptor
  except: return ""

def get_type(code):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  desc = cif.get("_chem_comp", {})[0]
  return getattr(desc, "type", "")

def get_name(code):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  desc = cif.get("_pdbx_chem_comp_identifier", {})
  for item in desc:
    if item.identifier:
      return item.identifier

def get_atom_names(code, alternate=False):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  tmp = []
  for item in cif["_chem_comp_atom"]:
    if alternate:
      tmp.append(item.alt_atom_id)
    else:
      tmp.append(item.atom_id)
  return tmp

def get_hydrogen_names(code,
                       wrap=False,
                       alternate=False
                       ):
  if wrap and alternate: assert 0
  cif = get_cif_dictionary(code)
  if not cif: return cif
  tmp = []
  for item in cif["_chem_comp_atom"]:
    if item.type_symbol in ["H"]:
      if wrap:
        name = item.atom_id
        if len(name)==4:
          name = name[3]+name[:3]
        tmp.append(name)
      elif alternate:
        tmp.append(item.alt_atom_id)
      else:
        tmp.append(item.atom_id)
  return tmp

def get_bond_pairs(code, alternate=False):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  tmp = []
  bonds = cif.get("_chem_comp_bond", {})
  for item in bonds:
    if alternate:
      atom1 = get_alternative_name(code, item.atom_id_1)
      atom2 = get_alternative_name(code, item.atom_id_2)
      tmp.append([atom1, atom2])
    else:
      tmp.append([item.atom_id_1, item.atom_id_2])
  return tmp

if __name__=="__main__":
  print '\nSMILES'
  print get_smiles(sys.argv[1])
  print '\nName'
  print get_name(sys.argv[1])
  print '\nAtom names'
  print get_atom_names(sys.argv[1])
  print '\nAlternate atom names'
  print get_atom_names(sys.argv[1], alternate=True)
  print '\nHydrogen names'
  print get_hydrogen_names(sys.argv[1])
  print '\nAlternate hydrogen names'
  print get_hydrogen_names(sys.argv[1], alternate=True)
  print '\nWrapped hydrogen names'
  print get_hydrogen_names(sys.argv[1], wrap=True)
  print '\nBond pairs'
  print get_bond_pairs(sys.argv[1])
  print '\nAlternate name bond pairs'
  print get_bond_pairs(sys.argv[1], alternate=True)
