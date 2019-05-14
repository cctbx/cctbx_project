
from __future__ import absolute_import, division, print_function
from libtbx.str_utils import make_sub_header
from libtbx.utils import Sorry, null_out
from libtbx import group_args, adopt_init_args, Auto, slots_getstate_setstate,\
      slots_getstate_setstate_default_initializer
import libtbx.phil
from math import sqrt
import sys

#-----------------------------------------------------------------------
# MODEL UTILITIES

def multi_conformer_selection(pdb_hierarchy):
  from scitbx.array_family import flex
  atoms = pdb_hierarchy.atoms()
  selection = flex.size_t()
  assert (not atoms.extract_i_seq().all_eq(0))
  for chain in pdb_hierarchy.only_model().chains( ):
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      if ((len(atom_groups) > 1) or (atom_groups[0].altloc.strip() != '')):
        selection.extend(residue_group.atoms().extract_i_seq())
  return selection

def fragment_single_conformer_chain(residues, chain_break_distance=3.0):
  """
  Split a protein chain into continuous peptide fragments (as lists of
  residue_group objects).
  """
  from scitbx.matrix import col
  k = 0
  current_fragment = []
  fragments = [ current_fragment ]
  while (k < len(residues)):
    res = residues[k]
    resseq = res.resseq_as_int()
    next_res = prev_res = None
    if (k > 0):
      prev_res = residues[k-1]
    if (k < len(residues) - 1):
      next_res = residues[k+1]
    k += 1
    if (prev_res is not None):
      prev_resseq = prev_res.resseq_as_int()
      c_atom = n_atom = None
      for atom in prev_res.atom_groups()[0].atoms():
        if (atom.name.strip() == "C"):
          c_atom = col(atom.xyz)
          break
      for atom in res.atom_groups()[0].atoms():
        if (atom.name.strip() == "N"):
          n_atom = col(atom.xyz)
          break
      if ((resseq > prev_resseq + 1) or
          (None in [n_atom, c_atom]) or
          (abs(n_atom-c_atom) > chain_break_distance)):
        current_fragment = [ res ]
        fragments.append(current_fragment)
        continue
    current_fragment.append(res)
  return fragments

def get_selection_gap(sel1, sel2):
  """
  Compute the gap or overlap between two selections (order-independent).
  Returns 0 if the selections are directly adjacent, or the number of
  residues overlapped (as a negative number) or missing between the
  selections (positive).
  """
  if (type(sel1).__name__ == 'bool'):
    sel1 = sel1.iselection()
  if (type(sel2).__name__ == 'bool'):
    sel2 = sel2.iselection()
  if (sel1[-1] == sel2[0] -1) or (sel2[-1] == sel1[0] - 1):
    return 0
  elif (sel2[-1] >= sel1[0] >= sel2[0]):
    return sel1[0] - sel2[-1] - 1
  elif (sel1[-1] >= sel2[0] >= sel1[0]):
    return sel2[0] - sel1[-1] - 1
  else :
    return max(sel2[0] - sel1[-1] - 1, sel1[0] - sel2[-1] - 1)

def score_rotamers(hierarchy, selection):
  """
  Count the number of rotamer outliers from a selection of residues in a
  PDB hierarchy.
  """
  from mmtbx.rotamer import rotamer_eval
  r = rotamer_eval.RotamerEval(data_version="8000")
  n_outliers = 0
  sub_hierarchy = hierarchy.select(selection) # XXX probably inefficient
  for rg in sub_hierarchy.only_model().only_chain().residue_groups():
    rotamer_flag = r.evaluate_residue(rg.only_atom_group())
    if (rotamer_flag == "OUTLIER"):
      n_outliers += 1
  return n_outliers

symmetric_atom_names = [
  ("OD1", "OD2"),
  ("OE1", "OE2"),
  ("OD1", "ND2"),
  ("OE1", "NE2"),
  ("CD1", "CD2"),
  ("CE1", "CE2"),
]
symmetric_atom_names_dict = dict(symmetric_atom_names +
    [ (n2,n1) for (n1,n2) in symmetric_atom_names ])

def coord_stats_with_flips(sites1, sites2, atoms):
  """
  Calculate RMSD and maximum distance for a pair of coordinate arrays,
  taking into account the symmetric or pseudo-symmetric nature of many
  sidechains.
  """
  from scitbx.matrix import col
  assert (len(sites1) == len(sites2) == len(atoms) > 0)
  rmsd_no_flip = rmsd_flip = None
  mean_sq = max_deviation_no_flip = 0
  n_sites = 0
  for site1, site2, atom in zip(sites1, sites2, atoms):
    if (atom.element.strip() == "H") : continue
    distance = abs(col(site1) - col(site2))
    mean_sq += distance**2
    if (distance > max_deviation_no_flip):
      max_deviation_no_flip = distance
    n_sites += 1
  assert (n_sites > 0)
  rmsd_no_flip = sqrt(mean_sq/n_sites)
  # TODO add HIS?
  if (not atoms[0].parent().resname in ["ASP","GLU","ASN","GLN","PHE","TYR"]):
    return group_args(rmsd=rmsd_no_flip, max_dev=max_deviation_no_flip)
  mean_sq = max_deviation_flip = 0
  for site1, site2, atom in zip(sites1, sites2, atoms):
    if (atom.element.strip() == "H") : continue
    atom_name = atom.name.strip()
    labels = atom.fetch_labels()
    symmetric_name = symmetric_atom_names_dict.get(atom_name, None)
    if (symmetric_name is not None):
      for site1_flip, site2_flip, atom_flip in zip(sites1, sites2, atoms):
        labels_flip = atom_flip.fetch_labels()
        if ((labels_flip.resid() == labels.resid()) and
            (labels_flip.chain_id == labels.chain_id) and
            (atom_flip.name.strip() == symmetric_name)):
          distance = abs(col(site1) - col(site2_flip))
          mean_sq += distance**2
          if (distance > max_deviation_flip):
            max_deviation_flip = distance
          break
      else : # didn't find the symmetry atom,
        rmsd_flip = float(sys.maxsize)
        max_deviation_flip = float(sys.maxsize)
        break
    else :
      distance = abs(col(site1) - col(site2))
      mean_sq += distance**2
      if (distance > max_deviation_flip):
        max_deviation_flip = distance
  rmsd_flip = sqrt(mean_sq/n_sites)
  return group_args(
    rmsd=min(rmsd_no_flip, rmsd_flip),
    max_dev=min(max_deviation_no_flip, max_deviation_flip))

def coord_stats_for_atom_groups(residue1, residue2):
  from scitbx.array_family import flex
  sites1 = flex.vec3_double()
  sites2 = flex.vec3_double()
  atoms = []
  for atom1 in residue1.atoms():
    if (atom1.element.strip() in ["H","D"]) : continue
    found = False
    for atom2 in residue2.atoms():
      if (atom2.name == atom1.name):
        assert (not found)
        found = True
        atoms.append(atom1)
        sites1.append(atom1.xyz)
        sites2.append(atom2.xyz)
  return coord_stats_with_flips(sites1, sites2, atoms)

def max_distance_between_atom_conformers(atom_groups, atom_name):
  atoms = []
  for atom_group in atom_groups :
    for atom in atom_group.atoms():
      if (atom.name.strip() == atom_name.strip()):
        atoms.append(atom)
  dxyz_max = 0
  if (len(atoms) > 1):
    for i_atm, atom1 in enumerate(atoms):
      for j_atom, atom2 in enumerate(atoms[(i_atm+1):]):
        dxyz = atom1.distance(atom2)
        if (dxyz > dxyz_max):
          dxyz_max = dxyz
  return dxyz_max

def requires_nterm_split(atom_groups, max_acceptable_distance=0.1):
  distance = max_distance_between_atom_conformers(atom_groups, "N")
  return distance > max_acceptable_distance

def requires_cterm_split(atom_groups, max_acceptable_distance=0.1):
  distance = max_distance_between_atom_conformers(atom_groups, "C")
  return distance > max_acceptable_distance

SEGID_NEW_REBUILT = "NEW1" # residues that have been fit to map
SEGID_NEW_SPLIT = "NEW2"   # residues that have been cloned without fitting
SEGID_MAIN = "OLD"         # original conformation of disordered residues
SELECTION_OLD = "segid OLD"
SELECTION_NEW_REBUILT = "segid NEW1"
SELECTION_NEW_SPLIT = "segid NEW2"
SELECTION_NEW = "segid NEW1 or segid NEW2"
SELECTION_MODIFIED = "segid NEW1 or segid NEW2 or segid OLD"

def spread_alternates(
    pdb_hierarchy,
    new_occupancy=None,
    split_all_adjacent=False,
    selection=None,
    verbose=False,
    log=None):
  """
  Given a model with some residues in alternate conformations, split the
  adjacent residues to allow for backbone movement.
  """
  if (log is None) : log = null_out()
  print("", file=log)
  print("  Splitting adjacent residues...", file=log)
  pdb_atoms = pdb_hierarchy.atoms()
  pdb_atoms.reset_i_seq()
  pdb_atoms.reset_tmp(0)
  if (selection is None):
    selection = pdb_hierarchy.atom_selection_cache().selection(
      SELECTION_NEW_REBUILT)
  elif isinstance(selection, str):
    selection = pdb_hierarchy.atom_selection_cache().selection(selection)
  def split_residue(residue_group):
    if (verbose):
      print("  %s" % residue_group.id_str(), file=log)
    new_occ = 0.5
    if (new_occupancy is not None):
      new_occ = max(0.2, min(0.8, new_occupancy))
    main_conf = residue_group.only_atom_group()
    main_conf.altloc = 'A'
    for atom in main_conf.atoms():
      atom.occ = 1.0 - new_occ
      atom.segid = SEGID_MAIN
    alt_conf = main_conf.detached_copy()
    alt_conf.altloc = 'B'
    for atom in alt_conf.atoms():
      atom.occ = new_occ
      atom.segid = SEGID_NEW_SPLIT
      atom.tmp = 1 # this flags the disordered segment for RSR
    residue_group.append_atom_group(alt_conf)
  n_split = 0
  for chain in pdb_hierarchy.only_model().chains():
    residue_groups = chain.residue_groups()
    for i_res, residue_group in enumerate(residue_groups):
      atom_groups = residue_group.atom_groups()
      resseq = residue_group.resseq_as_int()
      if (len(atom_groups) > 1):
        last_group_i_seqs = atom_groups[-1].atoms().extract_i_seq()
        if (not selection.select(last_group_i_seqs).all_eq(True)):
          #print "skipping %s" % residue_group.id_str()
          continue
        if (i_res > 0):
          prev_group = residue_groups[i_res - 1]
          prev_resseq = prev_group.resseq_as_int()
          if (prev_resseq >= resseq-1) and (len(prev_group.atom_groups())==1):
            if (split_all_adjacent) or (requires_nterm_split(atom_groups)):
              split_residue(prev_group)
              n_split += 1
        if (i_res < len(residue_groups) - 1):
          next_group = residue_groups[i_res + 1]
          next_resseq = next_group.resseq_as_int()
          if (next_resseq <= resseq+1) and (len(next_group.atom_groups())==1):
            if (split_all_adjacent) or (requires_cterm_split(atom_groups)):
              split_residue(next_group)
              n_split += 1
  return n_split

#-----------------------------------------------------------------------
# MAP STUFF
def get_partial_omit_map(
      fmodel,
      selection,
      selection_delete=None,
      negate_surrounding=False,
      map_file_name=None,
      partial_occupancy=0.5,
      resolution_factor=1/4.):
  """
  Generate an mFo-DFc map with a selection of atoms at reduced occupancy.
  Will write the map coefficients (along with 2mFo-DFc map) to an MTZ file
  if desired.  Reflections flagged for calculating R-free will always be
  omitted.
  """
  assert (0 <= partial_occupancy <= 1.0)
  xrs = fmodel.xray_structure
  occ = xrs.scatterers().extract_occupancies()
  occ.set_selected(selection, partial_occupancy)
  xrs.set_occupancies(occ)
  xrs_tmp = xrs.deep_copy_scatterers()
  if (selection_delete is not None):
    xrs_tmp = xrs_tmp.select(~selection_delete)
  fmodel.update_xray_structure(xrs_tmp, update_f_calc=True)
  fofc_coeffs = fmodel.map_coefficients(map_type="mFo-DFc",
    exclude_free_r_reflections=True)
  fofc_fft = fofc_coeffs.fft_map(resolution_factor=resolution_factor)
  fofc_map = fofc_fft.apply_sigma_scaling().real_map_unpadded()
  two_fofc_coeffs = fmodel.map_coefficients(map_type="2mFo-DFc",
    exclude_free_r_reflections=True)
  two_fofc_fft = two_fofc_coeffs.fft_map(resolution_factor=resolution_factor)
  two_fofc_map = two_fofc_fft.apply_sigma_scaling().real_map_unpadded()
  if (map_file_name is not None):
    import iotbx.map_tools
    iotbx.map_tools.write_map_coeffs(two_fofc_coeffs, fofc_coeffs,
      map_file_name)
  if (negate_surrounding):
    two_fofc_map = negate_surrounding_sites(
      map_data=two_fofc_map,
      xray_structure=xrs_tmp,
      iselection=selection)
    fofc_map = negate_surrounding_sites(
      map_data=fofc_map,
      xray_structure=xrs_tmp,
      iselection=selection)
  fmodel.update_xray_structure(xrs, update_f_calc=True)
  # XXX should the occupancies be reset now?
  return two_fofc_map, fofc_map

def negate_surrounding_sites(map_data, xray_structure, iselection,
      radius=1.5):
  """
  Makes the target map negative around an atom selection.
  """
  import mmtbx.refinement.real_space
  negate_selection = mmtbx.refinement.real_space.selection_around_to_negate(
    xray_structure          = xray_structure,
    selection_within_radius = 5, # XXX make residue dependent !!!!
    iselection              = iselection)
  target_map = mmtbx.refinement.real_space.\
    negate_map_around_selected_atoms_except_selected_atoms(
      xray_structure          = xray_structure,
      map_data                = map_data,
      negate_selection        = negate_selection,
      atom_radius             = radius)
  return target_map

filter_params_str = """
  discard_outliers = *rama *rota *cbeta *geo *map
    .type = choice(multi=True)
  min_model_map_cc = 0.85
    .type = float
  use_difference_map = True
    .type = bool
  sampling_radius = 2.5
    .type = float
  ignore_stub_residues = False
    .type = bool
"""

def is_validation_outlier(validation, params):
  filters = params.discard_outliers
  outlier = False
  if (validation.is_rotamer_outlier()) and ("rota" in filters):
    outlier = True
  if (validation.is_ramachandran_outlier()) and ("rama" in filters):
    outlier = True
  if (validation.is_cbeta_outlier()) and ("cbeta" in filters):
    outlier = True
  if (validation.is_geometry_outlier()) and ("geo" in filters):
    outlier = True
  if (validation.is_clash_outlier()) and ("clash" in filters):
    outlier = True
  if ((validation.is_map_outlier(cc_min=params.min_model_map_cc)) and
      ("map" in filters)):
    outlier = True
  return outlier

def filter_before_build(
    pdb_hierarchy,
    fmodel,
    geometry_restraints_manager,
    selection=None,
    params=None,
    verbose=True,
    log=sys.stdout):
  """
  Pick residues suitable for building alternate conformations - by default,
  this means no MolProbity/geometry outliers, good fit to map, no missing
  atoms, and no pre-existing alternates, but with significant difference
  density nearby.
  """
  from mmtbx.validation import molprobity
  from mmtbx.rotamer import rotamer_eval
  import mmtbx.monomer_library.server
  from mmtbx import building
  from iotbx.pdb import common_residue_names_get_class
  from scitbx.array_family import flex
  if (selection is None):
    selection = flex.bool(fmodel.xray_structure.scatterers().size(), True)
  pdb_atoms = pdb_hierarchy.atoms()
  assert (pdb_atoms.size() == fmodel.xray_structure.scatterers().size())
  pdb_atoms.reset_i_seq()
  full_validation = molprobity.molprobity(
    pdb_hierarchy=pdb_hierarchy,
    fmodel=fmodel,
    geometry_restraints_manager=geometry_restraints_manager,
    outliers_only=False,
    rotamer_library="8000")
  if (verbose):
    full_validation.show(out=log)
  multi_criterion = full_validation.as_multi_criterion_view()
  if (params is None):
    params = libtbx.phil.parse(filter_params_str).extract()
  mon_lib_srv = mmtbx.monomer_library.server.server()
  two_fofc_map, fofc_map = building.get_difference_maps(fmodel=fmodel)
  residues = []
  filters = params.discard_outliers
  make_sub_header("Identifying candidates for building", out=log)
  # TODO parallelize
  for chain in pdb_hierarchy.only_model().chains():
    if (not chain.is_protein()):
      continue
    for residue_group in chain.residue_groups():
      atom_groups = residue_group.atom_groups()
      id_str = residue_group.id_str()
      i_seqs = residue_group.atoms().extract_i_seq()
      residue_sel = selection.select(i_seqs)
      if (not residue_sel.all_eq(True)):
        continue
      if (len(atom_groups) > 1):
        print("  %s is already multi-conformer" % id_str, file=log)
        continue
      atom_group = atom_groups[0]
      res_class = common_residue_names_get_class(atom_group.resname)
      if (res_class != "common_amino_acid"):
        print("  %s: non-standard residue" % id_str, file=log)
        continue
      missing_atoms = rotamer_eval.eval_residue_completeness(
        residue=atom_group,
        mon_lib_srv=mon_lib_srv,
        ignore_hydrogens=True)
      if (len(missing_atoms) > 0):
        # residues modeled as pseudo-ALA are allowed by default; partially
        # missing sidechains are more problematic
        if ((building.is_stub_residue(atom_group)) and
            (not params.ignore_stub_residues)):
          pass
        else :
          print("  %s: missing or incomplete sidechain" % \
            (id_str, len(missing_atoms)), file=log)
          continue
      validation = multi_criterion.get_residue_group_data(residue_group)
      is_outlier = is_validation_outlier(validation, params)
      if (is_outlier):
        print("  %s" % str(validation), file=log)
        continue
      if (params.use_difference_map):
        i_seqs_no_hd = building.get_non_hydrogen_atom_indices(residue_group)
        map_stats = building.local_density_quality(
          fofc_map=fofc_map,
          two_fofc_map=two_fofc_map,
          atom_selection=i_seqs_no_hd,
          xray_structure=fmodel.xray_structure,
          radius=params.sampling_radius)
        if ((map_stats.number_of_atoms_below_fofc_map_level() == 0) and
            (map_stats.fraction_of_nearby_grid_points_above_cutoff()==0)):
          if (verbose):
            print("  no difference density for %s" % id_str, file=log)
          continue
      residues.append(residue_group.only_atom_group())
  if (len(residues) == 0):
    raise Sorry("No residues passed the filtering criteria.")
  print("", file=log)
  print("Alternate conformations will be tried for %d residue(s):" % \
      len(residues), file=log)
  building.show_chain_resseq_ranges(residues, out=log, prefix="  ")
  print("", file=log)
  return residues

rejoin_phil = libtbx.phil.parse("""
min_occupancy = 0.1
  .type = float
min_distance_over_u = 1.0
  .type = float
""")

class disorder_info(slots_getstate_setstate):
  __slots__ = [
    "n_confs", "sites", "atom_names", "id_str", "u_isos", "u_stars",
    "max_rmsd", "max_distance_over_u", "dxyz_save", "u_iso_save",
    "atom_name_save",
  ]
  def __init__(self, residue_group, xray_structure, ignore_hydrogens=True):
    from cctbx import adptbx
    from scitbx.array_family import flex
    from scitbx.matrix import col
    self.sites = []
    self.u_isos = []
    self.u_stars = []
    self.atom_names = []
    self.id_str = residue_group.id_str()
    unit_cell = xray_structure.unit_cell()
    atom_name_dict = {}
    for atom_group in residue_group.atom_groups():
      b_isos = atom_group.atoms().extract_b()
      if (len(b_isos) > 1) and (b_isos.all_eq(b_isos[0])):
        raise RuntimeError("B-factors for atom_group '%s' are identical" %
          atom_group.id_str())
      for atom in atom_group.atoms():
        if (atom.element.strip() in ["H", "D"]) and (ignore_hydrogens):
          continue
        if (not atom.name in atom_name_dict):
          atom_name_dict[atom.name] = []
        atom_name_dict[atom.name].append(atom)
    self.n_confs = 1
    for atom_name in sorted(atom_name_dict.keys()):
      n_atom_confs = len(atom_name_dict[atom_name])
      if (n_atom_confs > 1):
        if (self.n_confs > 1):
          if (self.n_confs != n_atom_confs):
            raise RuntimeError(("Inconsistent conformers for '%s': atom %s "+
              "has %d confs but previously %d conformations were seen") %
              (self.id_str, atom_name, n_atom_confs, self.n_confs))
        else :
          self.n_confs = n_atom_confs
          for n in range(n_atom_confs):
            self.sites.append(flex.vec3_double())
            self.u_isos.append(flex.double())
            self.u_stars.append([])
        self.atom_names.append(atom_name)
        for i_conf in range(n_atom_confs):
          atom = atom_name_dict[atom_name][i_conf]
          scatterer = xray_structure.scatterers()[atom.i_seq]
          assert scatterer.label == atom.id_str()
          self.sites[i_conf].append(atom.xyz)
          if (scatterer.flags.use_u_aniso()):
            self.u_isos[i_conf].append(adptbx.b_as_u(atom.b))
            self.u_stars[i_conf].append(scatterer.u_star)
          else :
            self.u_isos[i_conf].append(scatterer.u_iso)
            self.u_stars[i_conf].append(
              adptbx.u_iso_as_u_star(unit_cell, scatterer.u_iso))
    self.max_rmsd = 0
    self.max_distance_over_u = 0
    self.dxyz_save = self.u_iso_save = self.atom_name_save = None
    for i_conf, sites in enumerate(self.sites):
      if (i_conf < self.n_confs - 1):
        for j_conf in range(i_conf+1, self.n_confs):
          sites_next = self.sites[j_conf]
          rmsd = sites.rms_difference(sites_next)
          if (rmsd > self.max_rmsd):
            self.max_rmsd = rmsd
          for k_site, (site1, site2) in enumerate(zip(sites, sites_next)):
            distance = abs(col(site1) - col(site2))
            #u = max(self.u_isos[i_conf][k_site], self.u_isos[j_conf][k_site])
            proj_sum = adptbx.projection_sum(
              ustar1=self.u_stars[i_conf][k_site],
              ustar2=self.u_stars[j_conf][k_site],
              site1=site1,
              site2=site2,
              unit_cell=unit_cell)
            u_proj = max(proj_sum.z_12(), proj_sum.z_21())
            if (u_proj <= 0):
              continue
            ratio = distance / u_proj
            if (ratio > self.max_distance_over_u):
              self.max_distance_over_u = ratio
              self.dxyz_save = distance
              self.u_iso_save = u_proj
              self.atom_name_save = self.atom_names[k_site]

  def show_distance_info(self, out=sys.stdout, prefix="", suffix=''):
    if (suffix != '') : suffix = " " + suffix
    if (self.atom_name_save is not None):
      print(prefix + \
        "%s: rmsd=%-5.3f  atom '%s': dxyz=%-6.3f u=%-6.3f%s" % \
        (self.id_str, self.max_rmsd, self.atom_name_save, self.dxyz_save,
         sqrt(self.u_iso_save), suffix), file=out)
    else :
      print(prefix + "%s: rmsd=%-5.3f" % (self.id_str, self.max_rmsd), file=out)

def rejoin_split_single_conformers(
    pdb_hierarchy,
    params,
    crystal_symmetry=None,
    model_error_ml=None,
    reset_occupancies=Auto,
    verbose=False,
    log=sys.stdout):
  """
  Identifies residues which have been erroneously given alternate confomations
  which are actually identical - this is an expected side effect of the
  assembly process, but we want to keep the number of new atoms as limited as
  possible.  The filtering criteria are 1) the overall RMSD of all conformers
  of each residue to each other, and 2) the maximum ratio (over all non-H
  atoms) of the distance between conformations over u_iso (which in theory
  represents displacement in Angstroms).
  """
  from cctbx import adptbx
  from scitbx.array_family import flex
  from scitbx.matrix import col
  pdb_hierarchy.remove_hd()
  pdb_hierarchy.atoms().reset_i_seq()
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=crystal_symmetry)
  if (model_error_ml is None):
    model_error_ml = 0 # XXX ???
  n_modified = 0
  for chain in pdb_hierarchy.only_model().chains():
    residue_groups = chain.residue_groups()
    n_res = len(residue_groups)
    last_nconfs = 1
    pending_deletions = []
    for i_rg, residue_group in enumerate(residue_groups):
      residue_id = residue_group.id_str()
      atom_groups = residue_group.atom_groups()
      n_confs = len(atom_groups)
      if (n_confs > 1):
        info = disorder_info(residue_group=residue_group,
          xray_structure=xray_structure)
        delete_groups = set([])
        for i_group, atom_group in enumerate(atom_groups):
          occ = []
          for atom in atom_group.atoms():
            if (not atom.element.strip() in ["H","D"]):
              occ.append(atom.occ)
          if (len(occ) == 0) : # hydrogen-only
            continue
          if (max(occ) < params.min_occupancy):
            print("  %s: altloc '%s' occupancy = %4.2f" % (residue_id,
              atom_group.altloc, max(occ)), file=log)
            delete_groups.add(i_group)
        # TODO need to also check the C and N atoms specifically, because of
        # possible splitting of adjacent residues
        if ((info.max_rmsd < model_error_ml) or
            (info.max_distance_over_u < params.min_distance_over_u)):
          info.show_distance_info(out=log, prefix="  ", suffix='!!!')
          delete_groups.update(set(range(1, n_confs)))
        elif (verbose):
          info.show_distance_info(out=log, prefix="  ")
        reset_residue = False
        if (len(delete_groups) > 0):
          pending_deletions.append((i_rg, delete_groups))
    # now go back over the list of deletions and filter for continuous
    # fragments
    for k, (i_rg, delete_groups) in enumerate(pending_deletions):
      residue_group = residue_groups[i_rg]
      atom_groups = residue_group.atom_groups()
      n_confs = len(atom_groups)
      if (len(delete_groups) < n_confs - 1):
        kk = k + 1
        i_rg_last = i_rg
        while (kk < len(pending_deletions)):
          j_rg, delete_others = pending_deletions[kk]
          if (j_rg == i_rg_last + 1):
            delete_groups = set(delete_others).intersection(delete_groups)
          else :
            break
          kk += 1
      if (len(delete_groups) > 0):
        n_modified += 1
        #print >> log, "  merging residue %s" % atom_groups[0].id_str()
        for i_group, atom_group in enumerate(atom_groups):
          if (i_group in delete_groups):
            #print >> log, "    conformer %s deleted" % atom_group.altloc
            residue_group.remove_atom_group(atom_group)
            reset_residue = True
      # clumsy, but we're going to refine again anyway
      if (reset_residue):
        atom_groups = residue_group.atom_groups()
        for atom_group in atom_groups :
          for atom in atom_group.atoms():
            atom.occ = 1.0 / len(atom_groups)
          if (len(atom_groups) == 1):
            atom_group.altloc = ''
  if ((n_modified > 0) and reset_occupancies) or (reset_occupancies == True):
    atoms = pdb_hierarchy.atoms()
    occ = atoms.extract_occ()
    partial_occupancy_selection = occ < 1.0
    occ.set_selected(partial_occupancy_selection, 0.5)
    atoms.set_occ(occ)
  return n_modified

finalize_phil_str = """
set_b_iso = Auto
  .type = float
convert_to_isotropic = True
  .type = bool
"""

def finalize_model(pdb_hierarchy,
    xray_structure,
    set_b_iso=None,
    convert_to_isotropic=None,
    selection=None):
  """
  Prepare a rebuilt model for refinement, optionally including B-factor reset.
  """
  from cctbx import adptbx
  from scitbx.array_family import flex
  pdb_atoms = pdb_hierarchy.atoms()
  if (selection is None):
    selection = flex.bool(pdb_atoms.size(), True)
  elif isinstance(selection, str):
    sel_cache = pdb_hierarchy.atom_selection_cache()
    selection = sel_cache.selection(selection)
  for i_seq, atom in enumerate(pdb_atoms):
    assert (atom.parent() is not None)
    atom.segid = ""
    sc = xray_structure.scatterers()[i_seq]
    sc.label = atom.id_str()
  if (convert_to_isotropic):
    xray_structure.convert_to_isotropic(selection=selection.iselection())
  if (set_b_iso is not None):
    if (set_b_iso is Auto):
      u_iso = xray_structure.extract_u_iso_or_u_equiv()
      set_b_iso = adptbx.u_as_b(flex.mean(u_iso)) / 2
    xray_structure.set_b_iso(value=set_b_iso, selection=selection)
  pdb_hierarchy.adopt_xray_structure(xray_structure)
  pdb_atoms.reset_serial()
  pdb_atoms.reset_i_seq()
  pdb_atoms.reset_tmp()

class refined_fragment(slots_getstate_setstate_default_initializer):
  __slots__ = ["label", "selection", "sites_cart", "rmsd"]
  def show(self, out=sys.stdout, prefix=""):
    print(prefix+"refined %s (%d atoms): rmsd=%.3f" % (self.label,
      len(self.selection), self.rmsd), file=out)

class rsr_fragments_base(object):
  """
  Template for refinement of a collection of non-overlapping model fragments,
  intended to be used in parallel loops.
  """
  def __init__(self, pdb_hierarchy, fmodel, processed_pdb_file, **kwds):
    self.pdb_hierarchy = pdb_hierarchy
    self.fmodel = fmodel
    self.processed_pdb_file = processed_pdb_file
    self.__dict__.update(kwds)
    self.pdb_atoms = pdb_hierarchy.atoms()

  def make_box(self, selection, target_map):
    from mmtbx import building
    box = building.box_build_refine_base(
      pdb_hierarchy=self.pdb_hierarchy,
      xray_structure=self.fmodel.xray_structure,
      processed_pdb_file=self.processed_pdb_file,
      target_map=target_map,
      selection=selection,
      d_min=self.fmodel.f_obs().d_min(),
      out=null_out())
    return box

def print_trial_header(out, prefix=""):
  header_1 = """%20s%16s  %8s  %8s  %8s"""%("","  mFo-DFc map ","","","max")
  header_2 = """%-12s %5s  %7s  %7s  %8s  %8s  %8s""" % ("residue", "trial",
    "min." ,"mean", "CC", "RMSD", "change")
  print(prefix + header_1, file=out)
  print(prefix + header_2, file=out)
  print(prefix + "-" * (len(header_2)), file=out)

class trial_result(object):
  """
  Container for the results of a fitting run.  The sites_cart attribute
  should correspond to the iselection being refined.
  """
  def __init__(self, sites_cart, min_fofc, mean_fofc, cc, rmsd, max_dev):
    adopt_init_args(self, locals())

  def show_summary(self, prefix="", out=sys.stdout):
    print(prefix + ("min. mFo-DFc: %7.2f" % self.min_fofc), file=out)
    print(prefix + ("mean mFo-DFc: %7.2f" % self.mean_fofc), file=out)
    print(prefix + ("map-model CC: %8.3f" % self.cc), file=out)
    print(prefix + ("rmsd:         %8.3f" % self.rmsd), file=out)
    print(prefix + ("max. change:  %8.3f A" % self.max_dev), file=out)
    return self

  def __str__(self):
    return "%7.2f  %7.2f  %8.3f  %8.3f  %8.3f" % \
      (self.min_fofc, self.mean_fofc, self.cc, self.rmsd, self.max_dev)
