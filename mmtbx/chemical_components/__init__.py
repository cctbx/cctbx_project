
from __future__ import absolute_import, division, print_function
from mmtbx.chemical_components import cif_parser
from libtbx.utils import Sorry
import libtbx.load_env
import os
import sys
from six.moves import zip

l_amino_types = ['L-PEPTIDE LINKING']
d_amino_types = ['D-PEPTIDE LINKING']
amino_types = ['PEPTIDE LINKING',
               'PEPTIDE-LIKE',
               ] + l_amino_types + d_amino_types
rna_types = [
  'RNA LINKING',
  'L-RNA LINKING',
  'RNA OH 3 PRIME TERMINUS',
  'RNA OH 5 PRIME TERMINUS',
  ]
dna_types = [
  'DNA LINKING',
  'L-DNA LINKING',
  'DNA OH 3 PRIME TERMINUS',
  'DNA OH 5 PRIME TERMINUS',
  ]
rna_dna_types = rna_types + dna_types

sugar_types = ["SACCHARIDE",
               "L-SACCHARIDE",
               "D-SACCHARIDE",
               'D-SACCHARIDE, ALPHA LINKING',
               'D-SACCHARIDE, BETA LINKING',
               'L-SACCHARIDE, ALPHA LINKING',
               'L-SACCHARIDE, BETA LINKING',
               'D-SACCHARIDE 1,4 AND 1,4 LINKING',
               'L-SACCHARIDE 1,4 AND 1,4 LINKING',
               ]
terminii = [
  'L-PEPTIDE NH3 AMINO TERMINUS',
  'L-PEPTIDE COOH CARBOXY TERMINUS',
  'D-PEPTIDE NH3 AMINO TERMINUS',
]
non_polymer = [
  "NON-POLYMER",
  ]
non_alpha_peptide = [
  'L-BETA-PEPTIDE, C-GAMMA LINKING', # IAS
  'D-BETA-PEPTIDE, C-GAMMA LINKING', # ACB
  'D-GAMMA-PEPTIDE, C-DELTA LINKING', # FGA
  'L-GAMMA-PEPTIDE, C-DELTA LINKING', # GGL
  ]

loaded_cifs = {}

def find_data_dir():
  for relative_path in [
        "chem_data/chemical_components",
        "chemical_components",
        "ext_ref_files/chemical_components"]:
    result = libtbx.env.find_in_repositories(relative_path=relative_path)
    if (result is not None): return result
  return None

data_dir = find_data_dir()

def get_cif_filename(code):
  if (data_dir is None): return ""
  if (not code): return ""
  code=code.strip()
  if (len(code) == 0):
    raise Sorry("Residue code is blank.")
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

def get_cif_dictionary(code,
                       filename=None,
                       old_reader=False,
                       ):
  if filename is not None:
    if old_reader:
      cif = cif_parser.run2(filename)
    else:
      cif = cif_parser.run(filename)
  elif code in loaded_cifs:
    cif = loaded_cifs[code]
  else:
    filename = get_cif_filename(code)
    cif=None
    if os.path.exists(filename):
      if old_reader:
        cif = cif_parser.run2(filename)
      else:
        cif = cif_parser.run(filename)
      loaded_cifs[code] = cif
    if cif:
      for loop_id in ['_pdbx_chem_comp_descriptor', '_chem_comp']:
        for loop in cif.get(loop_id,[]):
          for attr in loop.__dict__:
            item = getattr(loop, attr)
            if type(item)==type(''):
              item = item.replace('\n', '')
              assert item.find('\n')==-1, '%s : %s' % (attr, item)
              setattr(loop, attr, str(item).strip())
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
      if item.program.strip()=="CACTVS" and item.type.strip()=="SMILES_CANONICAL":
        return item.descriptor
    except Exception: pass
  for item in desc:
    try:
      if item.program=="CACTVS" and item.type=="SMILES":
        return item.descriptor
    except Exception: pass
  try: return item.descriptor
  except Exception: return ""

def get_field_simple(code, loop, field):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  desc = cif.get(loop, {})[0]
  return getattr(desc, field, "")

def get_type(code):
  return get_field_simple(code, "_chem_comp", "type")

def get_name(code):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  desc = cif.get("_pdbx_chem_comp_identifier", {})
  for item in desc:
    if getattr(item, "identifier", None) is not None:
      return item.identifier
    else:
      return "Unknown or not read"

def get_atom_type_symbol(code):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  tmp = []
  for item in cif.get("_chem_comp_atom", []):
    tmp.append(item.type_symbol)
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

def get_atom_names(code, alternate=False, heavy_atom_only=False):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  tmp = []
  for item in cif["_chem_comp_atom"]:
    if heavy_atom_only:
      if item.type_symbol in ['H', 'D']: continue
    if alternate:
      tmp.append(item.alt_atom_id)
    else:
      tmp.append(item.atom_id)
  return tmp

def get_bond_pairs(code, alternate=False, heavy_atom_only=False, use_tuple=False):
  cif = get_cif_dictionary(code)
  if not cif: return cif
  if heavy_atom_only:
    heavy_atom_only = get_atom_names(code, heavy_atom_only=True)
  tmp = []
  bonds = cif.get("_chem_comp_bond", {})
  for item in bonds:
    if heavy_atom_only:
      if item.atom_id_1 not in heavy_atom_only or item.atom_id_2 not in heavy_atom_only:
        continue
    if alternate:
      atom1 = get_alternative_name(code, item.atom_id_1)
      atom2 = get_alternative_name(code, item.atom_id_2)
      pair = [atom1, atom2]
    else:
      pair = [item.atom_id_1, item.atom_id_2]
    pair.sort()
    if use_tuple:
      tmp.append(tuple(pair))
    else:
      tmp.append(pair)
  return tmp

def generate_chemical_components_codes(sort_reverse_by_smiles=False,
                                       only_non_polyer=False,
                                       ):
  def _cmp_smiles_length(f1, f2):
    c1 = f1[5:-4]
    c2 = f2[5:-4]
    if not c1:
      return 1
    if not c2:
      return -1
    s1 = get_smiles(c1)
    s2 = get_smiles(c2)
    if len(s1)>len(s2):
      return -1
    else:
      return 1

  data_dir = find_data_dir()
  dirs = os.listdir(data_dir)
  dirs.sort()
  filenames = []
  for d in dirs:
    if not os.path.isdir(os.path.join(data_dir, d)): continue
    filenames += os.listdir(os.path.join(data_dir, d))
  if sort_reverse_by_smiles:
    filenames.sort(_cmp_smiles_length)
  else:
    filenames.sort()
  for filename in filenames:
    if filename.find("data_")!=0: continue
    code = filename[5:-4]
    if only_non_polyer and get_group(code)!='non-polymer': continue
    yield code

def get_header(code):
  filename=get_cif_filename(code)
  if not filename: return ""
  f=open(filename)
  lines=f.readlines()
  f.close()
  outl=""
  for i, line in enumerate(lines):
    outl+="  %s" % line
    if i<3: continue
    if line.find("_chem_comp.")==-1:
      break
  return outl

def get_group(code, split_rna_dna=False, split_l_d=False):
  t = get_type(code)
  t=t.replace('"','').upper()
  if t in sugar_types:
    assert not split_l_d
    return 'sugar'
  elif t in amino_types:
    if split_l_d:
      if t in l_amino_types:
        return 'L-peptide'
      elif t in d_amino_types:
        return 'D-peptide'
      else:
        return 'amino_acid'
    return 'amino_acid'
  elif t in non_alpha_peptide:
    assert not split_l_d
    return 'non-alpha peptide'
  elif t in terminii:
    assert not split_l_d
    return 'amino_acid_terminal'
  elif t in rna_dna_types:
    assert not split_rna_dna
    if split_rna_dna:
      if t in rna_types:
        return 'rna'
      elif t in dna_types:
        return dna
      else:
        assert 0
    return 'rna_dna'
  elif t in ['NON-POLYMER',
             'PEPTIDE-LIKE',
            ]:
    return t.lower()
  print(t)
  assert 0

def get_restraints_group(code, split_rna_dna=True, split_l_d=True):
  g = get_group(code, split_rna_dna, split_l_d)
  print(g)
  if g in ['L-peptide', 'D-peptide', 'peptide',
           # 'non-polymer',
          ]:
    return g
  return {'amino_acid' : 'peptide',
          'non-polymer': 'ligand',
          }[g]
  assert 0

if __name__=="__main__":
  print('\nSMILES')
  print(get_smiles(sys.argv[1]))
  print('\nType')
  print(get_type(sys.argv[1]))
  print('\nWeight')
  print(get_field_simple(sys.argv[1], "_chem_comp", "formula_weight"))
  print('\nName')
  print(get_name(sys.argv[1]))
  print('\nAtom names')
  print(get_atom_names(sys.argv[1]))
  print('\nAlternate atom names')
  print(get_atom_names(sys.argv[1], alternate=True))
  print('\nHydrogen names')
  print(get_hydrogen_names(sys.argv[1]))
  print('\nAlternate hydrogen names')
  print(get_hydrogen_names(sys.argv[1], alternate=True))
  print('\nWrapped hydrogen names')
  print(get_hydrogen_names(sys.argv[1], wrap=True))
  print('\nBond pairs')
  print(get_bond_pairs(sys.argv[1]))
  print('\nAlternate name bond pairs')
  print(get_bond_pairs(sys.argv[1], alternate=True))
