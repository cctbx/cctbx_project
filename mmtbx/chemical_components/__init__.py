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

l_sugar_types = ["L-SACCHARIDE",
                 'L-SACCHARIDE, ALPHA LINKING',
                 'L-SACCHARIDE, BETA LINKING',
                 'L-SACCHARIDE 1,4 AND 1,4 LINKING',
                 ]
d_sugar_types = ["D-SACCHARIDE",
                 'D-SACCHARIDE, ALPHA LINKING',
                 'D-SACCHARIDE, BETA LINKING',
                 'D-SACCHARIDE 1,4 AND 1,4 LINKING',
                 ]
sugar_types = ["SACCHARIDE",
              ] + l_sugar_types + d_sugar_types
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

def is_chemical_components_file(filename):
  try:
    cif_model = iotbx.cif.reader(file_path=file_name).model()
    for cif_block in cif_model.values():
      if "_atom_site" in cif_block:
        return True
  except Exception as e:
    return False

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
                                       only_non_polymer=False,
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
    # if len(filenames)>100: break
  if sort_reverse_by_smiles:
    import functools
    filenames = sorted(filenames, key=functools.cmp_to_key(_cmp_smiles_length))
  else:
    filenames.sort()
  for filename in filenames:
    if filename.find("data_")!=0: continue
    code = filename[5:-4]
    if only_non_polymer and get_group(code)!='non-polymer': continue
    yield code

def get_header(code):
  filename=get_cif_filename(code)
  if not filename: return ""
  if not os.path.exists(filename): return ''
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

def get_group(code, split_rna_dna=False, split_l_d=False, verbose=False):
  t = get_type(code)
  if verbose: print('get_group',code, t)
  if t is not None:
    t=t.replace('"','').upper()
  else:
    if code in ['AD', 'TD']:
      t='L-DNA LINKING'
  if t in sugar_types:
    if split_l_d:
      if t in l_sugar_types:
        return 'L-saccharide'
      elif t in d_sugar_types:
        return 'D-saccharide'
    return 'saccharide'
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
    # assert not split_l_d
    return 'amino_acid_terminal'
  elif t in rna_dna_types:
    # assert not split_rna_dna
    if split_rna_dna:
      if t in rna_types:
        return 'RNA'
      elif t in dna_types:
        return 'DNA'
      else:
        assert 0
    return 'rna_dna'
  elif t in ['NON-POLYMER',
             'PEPTIDE-LIKE',
            ]:
    return t.lower()
  elif t in ['OTHER']:
    return t.lower()
  print('get_type for "%s" "%s"' % (code, t))
  assert 0

def get_restraints_group(code, split_rna_dna=True, split_l_d=True):
  g = get_group(code, split_rna_dna, split_l_d)
  print(g)
  if g in ['L-peptide', 'D-peptide', 'peptide',
            'L-saccharide', 'D-saccharide', 'saccharide',
           # 'non-polymer',
           'RNA', 'DNA',
          ]:
    return g
  return {'amino_acid'          : 'peptide',
          'amino_acid_terminal' : 'peptide',
          'non-polymer'         : 'ligand',
          # 'saccharide' : 'pyranose',
          }[g]
  assert 0

def get_as_atom_group(code):
  import iotbx
  from mmtbx.ligands.hierarchy_utils import _new_atom
  cif = get_cif_dictionary(code)
  if not cif: return cif
  tmp = []
  ag = iotbx.pdb.hierarchy.atom_group()
  ag.resname=code
  '''comp_id : BB9 <class 'str'>
  atom_id : HG <class 'str'>
  alt_atom_id : HG <class 'str'>
  type_symbol : H <class 'str'>
  charge : 0 <class 'int'>
  pdbx_align : 1 <class 'int'>
  pdbx_aromatic_flag : N <class 'str'>
  pdbx_leaving_atom_flag : N <class 'str'>
  pdbx_stereo_config : N <class 'str'>
  model_Cartn_x : 14.295 <class 'float'>
  model_Cartn_y : -4.046 <class 'float'>
  model_Cartn_z : 26.134 <class 'float'>
  pdbx_model_Cartn_x_ideal : -3.161 <class 'float'>
  pdbx_model_Cartn_y_ideal : -1.326 <class 'float'>
  pdbx_model_Cartn_z_ideal : -0.0 <class 'float'>
  pdbx_component_atom_id : HG <class 'str'>
  pdbx_component_comp_id : BB9 <class 'str'>
  pdbx_ordinal : 12 <class 'int'>'''
# def _new_atom(name, element, xyz, occ, b, hetero, segid=' '*4):
  use_model=True
  for item in cif["_chem_comp_atom"]:
    xyz = (item.model_Cartn_x,
           item.model_Cartn_y,
           item.model_Cartn_z,
           )
    print(xyz)
    if '?' in xyz:
      xyz = (item.pdbx_model_Cartn_x_ideal,
             item.pdbx_model_Cartn_y_ideal,
             item.pdbx_model_Cartn_z_ideal,
             )
      print(xyz)
      use_model=False
      break
  for item in cif["_chem_comp_atom"]:
    if use_model:
      xyz = (item.model_Cartn_x,
             item.model_Cartn_y,
             item.model_Cartn_z,
             )
    else:
      xyz = (item.pdbx_model_Cartn_x_ideal,
             item.pdbx_model_Cartn_y_ideal,
             item.pdbx_model_Cartn_z_ideal,
             )
    assert '?' not in xyz
    atom = _new_atom(item.atom_id,
                     item.type_symbol,
                     xyz,
                     1.,
                     20.,
                     True,
                     )
    ag.append_atom(atom)
  return ag

def get_as_hierarchy(code):
  import iotbx
  ag = get_as_atom_group(code)
  rg = iotbx.pdb.hierarchy.residue_group()
  rg.resseq='1'
  rg.append_atom_group(ag)
  chain = iotbx.pdb.hierarchy.chain()
  chain.id='A'
  chain.append_residue_group(rg)
  model = iotbx.pdb.hierarchy.model()
  model.append_chain(chain)
  ph = iotbx.pdb.hierarchy.root()
  ph.append_model(model)
  ph.reset_atom_i_seqs()
  return ph

def get_obsolete_status_from_chemical_components(code):
  # _chem_comp.pdbx_release_status                   OBS
  # _chem_comp.pdbx_replaced_by                      PCA
  header = get_header(code)
  pdbx_release_status = None
  pdbx_replaced_by = None
  if header is None: return None, None
  for line in header.split("\n"):
    tmp = line.split()
    if line.find("pdbx_release_status")>-1:
      pdbx_release_status = tmp[1]
    if line.find("pdbx_replaced_by")>-1:
      pdbx_replaced_by = tmp[1]
  return pdbx_release_status, pdbx_replaced_by

def get_filenames_from_start(start):
  rc=[]
  for code in generate_chemical_components_codes():
    if code.upper().find(start.upper())==0:
      rc.append(code)
  return rc

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
