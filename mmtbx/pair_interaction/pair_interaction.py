from __future__ import division
import os
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

def run_prepared(prepared, core=None, grid=None):
  """
  run(ph, core) from prepare(ph).

  Parameters
  ----------
  prepared : libtbx.group_args
      The result of prepare(ph).
  core : list of int, optional
      1-based residue group indices of the core; None for the whole
      hierarchy.
  grid : libtbx.group_args, optional
      xyz_min, xyz_max (Bohr): the grid of the searches instead of the
      bounding box of each search's atoms (_grid_search); the grid several
      views share (run_prepared_views, views_grid).

  Returns
  -------
  list or tuple
      As run: the interacting residue pairs (core None), else (core atoms,
      fragment atoms, fragment residues), atoms as 1-based positions in ph.
  """
  ph = prepared.ph
  atom_in_residue = prepared.atom_in_residue
  element_wfc_dict = prepared.element_wfc_dict
  if(core is None):
    return get_interactions(ph, list(atom_in_residue),
      element_wfc_dict=element_wfc_dict, grid=grid)
  cs = _core_selection(prepared, core)
  core_atoms = cs.core_atoms
  sub_ph = cs.sub_ph
  atom_in_residue_selected = cs.atom_in_residue
  #
  interactions = get_interactions(sub_ph, atom_in_residue_selected,
    silva_type='sedd', core=core, element_wfc_dict=element_wfc_dict,
    grid=grid)
  new_core = []
  for pair in interactions:
    if(len(set(pair).intersection(set(core))) > 0):
      new_core += pair
  new_core = list(set(new_core) | set(core))
  interactions = get_interactions(sub_ph, atom_in_residue_selected,
    silva_type='dori', core=new_core, element_wfc_dict=element_wfc_dict,
    grid=grid)
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

def _core_selection(prepared, core):
  """
  The part of a prepared hierarchy a core is searched in: core_atoms (i_seqs
  of the core residues, flex.size_t, in residue order), selection (flex.bool
  over the atoms: within 15 Bohr of a core atom, the only atoms that can
  interact with it), sub_ph (the selected atoms) and atom_in_residue (list,
  of the selected atoms). Shared by run_prepared and run_prepared_views so
  that the two cannot diverge.
  """
  core_atoms = flex.size_t()
  for idx in sorted(set(core)):
    if(idx <= prepared.n_residues):
      core_atoms.extend(prepared.residue_atoms[idx - 1])
  xyz = prepared.xyz
  selection = flex.bool(xyz.size(), False)
  for site in xyz.select(core_atoms):
    selection |= (xyz - site).norms() < 15
  return group_args(
    core_atoms      = core_atoms,
    selection       = selection,
    sub_ph          = prepared.ph.select(selection),
    atom_in_residue = list(prepared.atom_in_residue.select(selection)))

def _touching(pairs, residues):
  """The pairs (tuples) with a residue in residues (a set)."""
  return set([p for p in pairs if p[0] in residues or p[1] in residues])

def _view_cores(views, core):
  """
  Of run_prepared_views(views, core): to_view (per view, parent residue ->
  view residue), core_v (per view, the core in view residue indices),
  active (the views with a core residue) and sel (per active view, its
  _core_selection).
  """
  to_view, core_v = [], []
  for v in views:
    d = dict([(r, k+1) for k, r in enumerate(v.residues)])
    to_view.append(d)
    core_v.append([d[r] for r in core if r in d])
  active = [i for i in range(len(views)) if len(core_v[i]) > 0]
  sel = dict([(i, _core_selection(views[i].prepared, core_v[i]))
              for i in active])
  return group_args(to_view=to_view, core_v=core_v, active=active, sel=sel)

def _union_box(views, vc):
  """group_args(xyz_min, xyz_max): the box of the searched atoms of all
  active views (Bohr); None without an active view."""
  xyz = flex.vec3_double()
  for i in vc.active:
    xyz.extend(views[i].prepared.xyz.select(vc.sel[i].selection))
  if(xyz.size() == 0): return None
  return group_args(xyz_min=xyz.min(), xyz_max=xyz.max())

def views_grid(views, core=None):
  """
  The grid run_prepared_views(views, core) searches every view on: the
  bounding box (Bohr) of the union over the views of the atoms within 15
  Bohr of the core (_core_selection); without a core (views_pairs), the
  box of all atoms of all views. Views of duplicate conformers (same
  coordinates) get the box run_prepared alone would use, so their results
  are those of the single-conformer hierarchy, bit for bit.
  None when no view has a residue of the core (or an atom).
  """
  if(core is None):
    xyz = flex.vec3_double()
    for v in views: xyz.extend(v.prepared.xyz)
    if(xyz.size() == 0): return None
    return group_args(xyz_min=xyz.min(), xyz_max=xyz.max())
  return _union_box(views, _view_cores(views, core))

def run_prepared_views(views, core):
  """
  run_prepared of several single-conformer views of one parent hierarchy
  for a core given in parent residue indices, every view on one shared
  grid (views_grid): the first view with a core residue (the reference) is
  searched in full, the others only at the grid points within 10 Bohr of
  an atom that differs between the views; a far point sees the same atoms
  in every view and gives the same pair, so every view's result equals its
  independent run_prepared(grid=views_grid(views, core)).

  Parameters
  ----------
  views : list of libtbx.group_args
      prepared: prepare(view);
      iselection: flex.size_t, the parent atom position of each view atom
      (ascending: a view is a selection of the parent);
      residues: list, the 1-based parent residue group index of each view
      residue group.
  core : list of int
      1-based parent residue group indices.

  Returns
  -------
  list
      Per view: None when the view has no residue of the core, else
      (core atoms, fragment atoms, fragment residues), sorted lists in
      parent terms (1-based parent atom positions, parent residue indices),
      equal to run_prepared(view.prepared, core in view indices, grid)
      mapped to the parent.
  """
  n = len(views)
  vc = _view_cores(views, core)
  to_view, core_v, active, sel = vc.to_view, vc.core_v, vc.active, vc.sel
  result = [None]*n
  if(len(active) == 0): return result
  ref = active[0]
  grid = _union_box(views, vc)
  # the searched part of every view, as parent atoms
  parent_atoms = dict([(i, views[i].iselection.select(sel[i].selection))
                       for i in active])
  parent_sets = dict([(i, set(parent_atoms[i])) for i in active])
  # the differing atoms (region): symmetric differences to the reference
  differing = set()
  for i in active:
    if(i != ref): differing |= parent_sets[i] ^ parent_sets[ref]
  region = flex.vec3_double()
  if(len(differing) > 0):
    seen = set()
    for i in active:
      xyz = views[i].prepared.xyz.select(sel[i].selection)
      for k, a in enumerate(parent_atoms[i]):
        if(a in differing and a not in seen):
          seen.add(a)
          region.append(xyz[k])
  if(region.size() == 0): region = None
  core_p = dict([(i, set([views[i].residues[k-1] for k in core_v[i]]))
                 for i in active])

  def search(i, silva_type, core_i, near_only):
    """Pairs of view i in parent residue indices, with their flags."""
    step = {'sedd': 0.3, 'dori': 0.5}[silva_type]*A2B
    d = _grid_search(sel[i].sub_ph, sel[i].atom_in_residue, step,
      silva_type, core_i, views[i].prepared.element_wfc_dict, grid=grid,
      region_sites=region, near_only=near_only)
    res = views[i].residues
    out = {}
    for (a, b), flag in d.items():
      pa, pb = res[a-1], res[b-1]
      pair = (min(pa, pb), max(pa, pb))
      if(flag < out.get(pair, 2)): out[pair] = flag
    return out

  # pass 1 (sedd): the reference in full, the others at near points
  d_ref = search(ref, 'sedd', core_v[ref], False)
  far_ref = set([p for p, f in d_ref.items() if f == 0])
  new_core = {}
  for i in active:
    if(i == ref):
      pairs = set(d_ref)
    else:
      pairs = _touching(far_ref, core_p[i])
      if(region is not None):
        pairs |= set(search(i, 'sedd', core_v[i], True))
    nc = set(core_p[i])
    for p in _touching(pairs, core_p[i]): nc.update(p)
    new_core[i] = nc
  # pass 2 (dori): the reference with the union of the new cores
  union = set()
  for i in active: union |= new_core[i]
  core2_ref = sorted([to_view[ref][r] for r in union if r in to_view[ref]])
  d2_ref = search(ref, 'dori', core2_ref, False)
  far2_ref = set([p for p, f in d2_ref.items() if f == 0])
  for i in active:
    if(i == ref):
      pairs = _touching(set(d2_ref), new_core[i])
    else:
      pairs = _touching(far2_ref, new_core[i])
      if(region is not None):
        core2_i = sorted([to_view[i][r] for r in new_core[i]])
        pairs |= set(search(i, 'dori', core2_i, True))
    mols = set(new_core[i])
    for p in _touching(pairs, core_p[i]): mols.update(p)
    isel = views[i].iselection
    residue_atoms = views[i].prepared.residue_atoms
    atoms = []
    for m in mols:
      atoms.extend([isel[j]+1 for j in residue_atoms[to_view[i][m]-1]])
    result[i] = (sorted([isel[j]+1 for j in sel[i].core_atoms]),
                 sorted(atoms), sorted(mols))
  return result

def interaction_pairs(prepared, silva_type, grid=None):
  """
  The interacting residue pairs of a whole prepared hierarchy.

  Parameters
  ----------
  prepared : libtbx.group_args
      The result of prepare(ph).
  silva_type : str
      'sedd' (step 0.3 A: the pass of run_prepared that extends the core)
      or 'dori' (step 0.5 A: the pass that adds the buffer).
  grid : libtbx.group_args, optional
      xyz_min, xyz_max (Bohr), as run_prepared.

  Returns
  -------
  set of tuple
      (i, j), i < j, 1-based residue group indices in ph.
  """
  step = {'sedd': 0.3, 'dori': 0.5}[silva_type]*A2B
  return set(_grid_search(prepared.ph, list(prepared.atom_in_residue), step,
    silva_type, None, prepared.element_wfc_dict, grid=grid))

def partner_map(pairs):
  """{residue: set of its partners} of a set of pairs (interaction_pairs)."""
  result = {}
  for a, b in pairs:
    result.setdefault(a, set()).add(b)
    result.setdefault(b, set()).add(a)
  return result

def fragment_residues(core, sedd, dori):
  """
  The rule of run_prepared applied to the pairs of the whole hierarchy:
  new_core, the core and its sedd partners; fragment, new_core and the
  dori partners of the core. The points that give a pair touching a core
  residue lie within 10 Bohr of a core atom, so on one grid the
  core-restricted search of run_prepared finds the same pairs as the
  whole-hierarchy search. run_prepared alone differs through its 15 Bohr
  sub-selection around the core, which truncates the environment of the
  points near its edge and sets the grid origin (its bounding box). One
  search per pass for the whole hierarchy then serves every core
  (fragmentation: one super sphere, many clusters).

  Parameters
  ----------
  core : list of int
      1-based residue group indices.
  sedd, dori : dict
      partner_map of the sedd and dori pairs.

  Returns
  -------
  tuple of set
      (new_core, fragment): 1-based residue group indices.
  """
  new_core = set(core)
  for r in core: new_core |= sedd.get(r, set())
  fragment = set(new_core)
  for r in core: fragment |= dori.get(r, set())
  return new_core, fragment

def fragment_from_pairs(prepared, core, sedd_pairs, dori_pairs):
  """
  run_prepared(prepared, core) from the pairs of the whole hierarchy
  (interaction_pairs, 'sedd' and 'dori'; fragment_residues): (core atoms,
  fragment atoms, fragment residues), sorted lists, atoms 1-based positions
  in ph. A core residue beyond the hierarchy has no atoms.
  """
  new_core, fragment = fragment_residues(core, partner_map(sedd_pairs),
                                         partner_map(dori_pairs))
  def atoms(residues):
    result = []
    for r in residues:
      if(r <= prepared.n_residues):
        result.extend([i+1 for i in prepared.residue_atoms[r-1]])
    return sorted(result)
  return (atoms(set(core)), atoms(fragment), sorted(fragment))

def views_pairs(views, silva_type):
  """
  interaction_pairs of several single-conformer views of one parent
  hierarchy, united in parent residue indices, every view on one shared
  grid (views_grid(views)): the first view is searched in full, the others
  only at the grid points within 10 Bohr of an atom that differs between
  the views; a far point sees the same atoms in every view and gives the
  same pair (run_prepared_views, with the whole views as the searched
  part), so the result is the union of the views' independent searches on
  that grid.

  Parameters
  ----------
  views : list of libtbx.group_args
      prepared, iselection, residues, as run_prepared_views; a view
      without atoms contributes nothing.
  silva_type : str
      As interaction_pairs.

  Returns
  -------
  set of tuple
      (i, j), i < j, 1-based parent residue group indices.
  """
  result = set()
  views = [v for v in views if v.prepared.xyz.size() > 0]
  if(len(views) == 0): return result
  grid = views_grid(views)
  # the differing atoms (region): the parent atoms not in every view
  n_parent = max([flex.max(v.iselection) for v in views]) + 1
  count = flex.int(n_parent, 0)
  for v in views:
    count.set_selected(v.iselection, count.select(v.iselection) + 1)
  differing = count < len(views)
  region = flex.vec3_double()
  seen = flex.bool(n_parent, False)
  for v in views:
    isel = v.iselection
    take = differing.select(isel) & ~seen.select(isel)
    region.extend(v.prepared.xyz.select(take))
    seen.set_selected(isel.select(take), True)
  if(region.size() == 0): region = None
  step = {'sedd': 0.3, 'dori': 0.5}[silva_type]*A2B
  for i, v in enumerate(views):
    if(i > 0 and region is None): break
    d = _grid_search(v.prepared.ph, list(v.prepared.atom_in_residue), step,
      silva_type, None, v.prepared.element_wfc_dict, grid=grid,
      region_sites=region, near_only=(i > 0))
    res = v.residues
    for a, b in d:
      pa, pb = res[a-1], res[b-1]
      result.add((min(pa, pb), max(pa, pb)))
  return result

def get_interactions(ph, atom_in_residue, step_size=0.5*A2B, silva_type='dori',
      core=None, element_wfc_dict=None, grid=None):
  """
  The interacting residue pairs (1-based residue indices of atom_in_residue)
  of ph on the grid given (_grid_search), else on its own bounding box:
  _grid_search without flags.
  element_wfc_dict: the wave functions by element key (prepare).
  """
  if(silva_type=='sedd'):step_size=0.3*A2B
  pairs = _grid_search(ph, atom_in_residue, step_size, silva_type, core,
    element_wfc_dict, grid=grid)
  return list(set(pairs.keys()))

def _grid_search(ph, atom_in_residue, step_size, silva_type, core,
      element_wfc_dict, grid=None, region_sites=None, near_only=False):
  """
  ext.points_and_pairs of ph with the flags.

  Parameters
  ----------
  ph : iotbx.pdb.hierarchy.root
      Coordinates in Bohr (prepare).
  atom_in_residue : list of int
      1-based residue index of every atom of ph.
  step_size : float
      Grid spacing (Bohr).
  silva_type : str
      'sedd' or 'dori'.
  core : list of int or None
      Residue indices whose pairs are wanted; None for all pairs.
  element_wfc_dict : dict
      Wave functions by element key (prepare).
  grid : libtbx.group_args, optional
      xyz_min, xyz_max (Bohr): the grid origin and extent, shared between
      hierarchies; every atom must lie inside. Default: the bounding box
      of the atoms of ph.
  region_sites : flex.vec3_double, optional
      Coordinates (Bohr) of the atoms that differ between the hierarchies
      sharing the grid; a grid point within 10 Bohr of one is near.
  near_only : bool
      Search the near points only.

  Returns
  -------
  dict
      {(i, j): flag}, i < j 1-based residue indices; flag 0 when the pair
      was found at a far point (always, without region_sites), 1 when it
      was found only at near points.
  """
  atoms = ph.atoms()
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
  xyz = atoms.extract_xyz()
  if(grid is None):
    xyz_min = xyz.min()
    xyz_max = xyz.max()
  else:
    xyz_min = tuple(grid.xyz_min)
    xyz_max = tuple(grid.xyz_max)
    lo = xyz.min()
    hi = xyz.max()
    for k in range(3):
      assert lo[k] >= xyz_min[k] and hi[k] <= xyz_max[k], (
        "atom outside the shared grid: axis %d, atoms %s, grid %s" % (
          k, (lo[k], hi[k]), (xyz_min[k], xyz_max[k])))
  xyz_step = [
    int(math.floor((xyz_max[i]-xyz_min[i])/step_size)+1) for i in range(3)]
  # Only pairs with a core residue are used by run(): the others are skipped
  residue_in_core = flex.bool()
  if(core is not None):
    residue_in_core = flex.bool(max(atom_in_residue)+1, False)
    for i in core:
      if(i < residue_in_core.size()): residue_in_core[i] = True
  if(region_sites is None): region_sites = flex.vec3_double()
  interacting_pairs = ext.points_and_pairs(
    ngrid           = xyz_step,
    step_size       = step_size,
    xyz             = xyz,
    xyz_min         = xyz_min,
    atom_in_residue = atom_in_residue,
    element_flags   = element_flags,
    wfc_obj         = wave_functions,
    silva_type      = silva_type,
    residue_in_core = residue_in_core,
    region_sites    = region_sites,
    near_only       = near_only)
  result = {}
  for it in interacting_pairs:
    pair = (int(it[0]), int(it[1]))
    flag = int(it[2])
    if(flag < result.get(pair, 2)): result[pair] = flag
  return result
