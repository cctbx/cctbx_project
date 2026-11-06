from __future__ import division
import os
import time
import math
import libtbx.load_env
from scitbx.array_family import flex
from libtbx import easy_pickle
from libtbx.utils import Sorry
from libtbx import group_args

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("mmtbx_pair_interaction_ext")

dat_path = libtbx.env.find_in_repositories("qrefine")
qr_unit_tests_data = None
if(dat_path is not None):
  qr_unit_tests_data = os.path.join(dat_path,"tests","unit","data_files")

A2B=1.8897259885789

_WFC_CACHE = {} # key -> ext.wfc, per process

def _wfc_key(element):
  """'c_' for C, 'cl' for Cl: the element symbol, lower case, padded with _."""
  key = element.strip().lower()
  if(len(key) == 1): key += "_"
  return key

def load_wfc(element, folder=None):
  """
  Wave function object of an element.

  Parameters
  ----------
  element : str
      Element symbol in any case; the key form ('c_', 'cl') is accepted.
  folder : str, optional
      Folder of the pickles <key>_wfc_obj.pkl; mmtbx/pair_interaction by
      default.

  Returns
  -------
  mmtbx_pair_interaction_ext.wfc
      Loaded from <folder>/<key>_wfc_obj.pkl, else parsed from the text
      table qrefine/plugin/yoink/dat/<key>_lda.wfc and written as that
      pickle into folder (silently not written when folder is not
      writable). Objects are cached per process. Nothing is written into
      the current directory.
  """
  key = _wfc_key(element)
  wfc_obj = _WFC_CACHE.get(key)
  if(wfc_obj is not None): return wfc_obj
  if(folder is None):
    folder = libtbx.env.find_in_repositories("mmtbx/pair_interaction")
  pkl = os.path.join(folder, key + "_wfc_obj.pkl")
  if(os.path.isfile(pkl)):
    wfc_obj = easy_pickle.load(pkl)
  else:
    wfc_obj = _parse_wfc(key)
    try:
      easy_pickle.dump(pkl, wfc_obj)
    except (IOError, OSError):
      pass
  _WFC_CACHE[key] = wfc_obj
  return wfc_obj

def _parse_wfc(key):
  """ext.wfc from the text table qrefine/plugin/yoink/dat/<key>_lda.wfc."""
  folder = libtbx.env.find_in_repositories("qrefine/plugin/yoink/dat")
  if(folder is None):
    raise Sorry("No _lda.wfc files found.")
  lines=open(os.path.join(folder, key + "_lda.wfc")).readlines()
  num_orbitals=int(lines[0].strip().split()[0])
  occ_electrons_array = [[0]*num_orbitals]
  line_three = lines[2].strip().split()
  for i in range(len(line_three)):
    occ_electrons_array[0][i]=int(line_three[i])
  occ_electrons=occ_electrons_array[0]

  line_four = lines[3].strip().split()
  xmin      = float(line_four[0])
  zz        = float(line_four[1])
  dx        = float(line_four[2])
  ngrid     = int(line_four[3])

  r_array = flex.double()
  wfcin_array  = []
  for line in lines[4:]:
    line = line.split()
    r_array.append(float(line[0]))
    tmp = []
    for i in range(1,len(line)):
      tmp.append(float(line[i]))
    wfcin_array.append(tmp)
  assert len(wfcin_array) == ngrid, [len(wfcin_array) , ngrid]

  return ext.wfc(
    ngrid         = ngrid,
    zz            = zz,
    r_array       = r_array,
    wfcin_array   = wfcin_array,
    occ_electrons = occ_electrons)

def prepare(ph):
  """
  Per-hierarchy data of run_prepared, computed once per hierarchy.

  Parameters
  ----------
  ph : iotbx.pdb.hierarchy.root
      The hierarchy; it is not changed.

  Returns
  -------
  libtbx.group_args
      ph: a deep copy with coordinates in Bohr and i_seqs reset; atoms: its
      atoms; xyz: its coordinates; atom_in_residue (flex.int): the 1-based
      residue group index of every atom; residue_atoms: the i_seqs
      (flex.size_t) of every residue group; n_residues; element_wfc_dict:
      the wave functions of the elements present, by key (load_wfc).
      Prepared once, a hierarchy serves run_prepared for every core
      (fragmentation: one super sphere, many clusters).
  """
  ph = ph.deep_copy()
  ph.atoms().reset_i_seq()
  atom_in_residue = flex.int()
  residue_atoms = []
  element_types = set()
  for i_rg, rg in enumerate(ph.residue_groups()):
    atoms = rg.atoms()
    residue_atoms.append(atoms.extract_i_seq())
    atom_in_residue.extend(flex.int(atoms.size(), i_rg + 1))
    for e in atoms.extract_element():
      element_types.add(_wfc_key(e))
  assert atom_in_residue.size() == ph.atoms_size()
  xyz = ph.atoms().extract_xyz() * A2B
  ph.atoms().set_xyz(xyz)
  element_wfc_dict = {}
  for element in element_types:
    element_wfc_dict[element] = load_wfc(element)
  return group_args(
    ph               = ph,
    atoms            = ph.atoms(),
    xyz              = xyz,
    atom_in_residue  = atom_in_residue,
    residue_atoms    = residue_atoms,
    n_residues       = len(residue_atoms),
    element_wfc_dict = element_wfc_dict)

def run(ph, core=None):
  """
  core=None: the list of interacting residue pairs (1-based residue group
  indices in ph) of the whole hierarchy.
  core=[residue indices]: (core atoms, fragment atoms, fragment residues) of
  the core residues, the fragment being the core plus the residues it
  interacts with. Atoms are 1-based positions in ph.
  ph is not changed (prepare works on a copy).
  """
  return run_prepared(prepare(ph), core)

def run_prepared(prepared, core=None):
  """
  run(ph, core) from prepare(ph).

  Parameters
  ----------
  prepared : libtbx.group_args
      The result of prepare(ph).
  core : list of int, optional
      1-based residue group indices of the core; None for the whole
      hierarchy.

  Returns
  -------
  list or tuple
      As run: the interacting residue pairs (core None), else (core atoms,
      fragment atoms, fragment residues), atoms as 1-based positions in ph.
  """
  ph = prepared.ph
  xyz = prepared.xyz
  atom_in_residue = prepared.atom_in_residue
  element_wfc_dict = prepared.element_wfc_dict
  if(core is None):
    return get_interactions(ph, list(atom_in_residue),
      element_wfc_dict=element_wfc_dict)
  # core atoms in residue order (as the residue groups come)
  core_atoms = flex.size_t()
  for idx in sorted(set(core)):
    if(idx <= prepared.n_residues):
      core_atoms.extend(prepared.residue_atoms[idx - 1])
  # Only atoms within 15 Bohr of the core can interact with it: the
  # density is evaluated for the core and these atoms.
  selection = flex.bool(xyz.size(), False)
  for site in xyz.select(core_atoms):
    selection |= (xyz - site).norms() < 15
  sub_ph = ph.select(selection)
  atom_in_residue_selected = list(atom_in_residue.select(selection))
  #
  interactions = get_interactions(sub_ph, atom_in_residue_selected,
    silva_type='sedd', core=core, element_wfc_dict=element_wfc_dict)
  new_core = []
  for pair in interactions:
    if(len(set(pair).intersection(set(core))) > 0):
      new_core += pair
  new_core = list(set(new_core) | set(core))
  interactions = get_interactions(sub_ph, atom_in_residue_selected,
    silva_type='dori', core=new_core, element_wfc_dict=element_wfc_dict)
  interaction_mols = new_core
  for item in interactions:
    if(len(set(item).intersection(set(core))) > 0):
      interaction_mols = interaction_mols + list(item)
  interaction_atoms = []
  for mol_id in interaction_mols:
    interaction_atoms += [i + 1 for i in prepared.residue_atoms[mol_id - 1]]
  return (list(set([i + 1 for i in core_atoms])),
          list(set(interaction_atoms)),
          list(set(interaction_mols)))

def get_interactions(ph, atom_in_residue, step_size=0.5*A2B, silva_type='dori',
      core=None, element_wfc_dict=None):
  # element_wfc_dict: the wave functions by element key (prepare)
  if(silva_type=='sedd'):step_size=0.3*A2B
  t0=time.time()
  atoms = ph.atoms()
  #print("num atoms:",len(atoms))
  elements = atoms.extract_element()
  eldict = dict(enumerate(set(elements)))
  tmp = {}
  for k, v in zip(eldict.keys(), eldict.values()):
    tmp[v.strip()] = k
  eldict = tmp
  element_flags = flex.int()
  for e in elements:
    element_flags.append(eldict[e.strip()])
  # element_flags index into wave_functions, so build it for this call only
  wave_functions = []
  for e in eldict:
    k = e.lower()
    if(len(k)==1): k+="_"
    wave_functions.append(element_wfc_dict[k])
  xyz = ph.atoms().extract_xyz()
  xyz_min=xyz.min()
  xyz_max=xyz.max()
  xyz_step = [
    int(math.floor((xyz_max[i]-xyz_min[i])/step_size)+1) for i in range(3)]
  #print("xyz_min,xyz_max",xyz_min,xyz_max)
  #print("step_size", step_size)
  #print("xyz_step", xyz_step)
  #print("points to process:",xyz_step[0]*xyz_step[1]*xyz_step[2])
  # Only pairs with a core residue are used by run(): the others are skipped
  residue_in_core = flex.bool()
  if(core is not None):
    residue_in_core = flex.bool(max(atom_in_residue)+1, False)
    for i in core:
      if(i < residue_in_core.size()): residue_in_core[i] = True
  interacting_pairs = ext.points_and_pairs(
    ngrid           = xyz_step,
    step_size       = step_size,
    xyz             = xyz,
    xyz_min         = xyz_min,
    atom_in_residue = atom_in_residue,
    element_flags   = element_flags,
    wfc_obj         = wave_functions,
    silva_type      = silva_type,
    residue_in_core = residue_in_core)
  #print("interacting_pairs", list(set(list((interacting_pairs)))) )
  tmp = []
  for it in interacting_pairs:
    pair = [int(it[0]),int(it[1])]
    tmp.append(tuple(pair))
  interacting_pairs = list(set(tmp))
  #print(interacting_pairs)
  #print("Time ",(time.time()-t0))
  return interacting_pairs
