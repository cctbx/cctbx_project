
# TODO more tests?

"""
Prototype for building alternate conformations into difference density.
The actual method is a variant of one that Pavel suggested, in combination
with the procedure in mmtbx.building.extend_sidechains: first, the
backbone atoms for a residue and its neighbors are refined into the
target map using minimization and/or annealing, then the sidechain is
replaced (using an idealized copy) and its placement optimized by a grid
search that also allows for backbone flexibility.
"""

from __future__ import absolute_import, division, print_function
from mmtbx.building import extend_sidechains
from mmtbx.building import alternate_conformations as alt_confs
from mmtbx import building
from libtbx import adopt_init_args, Auto, slots_getstate_setstate
from libtbx.str_utils import make_header, make_sub_header, format_value
from libtbx.utils import null_out
from libtbx import easy_mp
from cStringIO import StringIO
import time
import sys
from six.moves import range

build_params_str = """
  expected_occupancy = None
    .type = float(value_min=0.0,value_max=1.0)
  omit_waters = False
    .type = bool
  window_size = 2
    .type = int
  backbone_sample_angle = 10
    .type = int
  anneal = False
    .type = bool
  annealing_temperature = 1000
    .type = int
  simple_chi1_sampling = True
    .type = bool
  rmsd_min = 0.5
    .type = float
  map_thresholds {
    omit_two_fofc_min_sc_mean = 0.8
      .type = float
    omit_two_fofc_min_mc = 1.0
      .type = float
    omit_fofc_min_sc_mean = 3.0
      .type = float
    two_fofc_min = 1.0
      .type = float
    fofc_min = 3.0
      .type = float
  }
"""

master_phil_str = """
residue_fitting {
  %s
  #delete_hydrogens = False
  #  .type = bool
}
prefilter {
  include scope mmtbx.building.alternate_conformations.filter_params_str
}
cleanup {
  rsr_after_build = True
    .type = bool
  rsr_fofc_map_target = True
    .type = bool
  rsr_max_cycles = 3
    .type = int
  include scope mmtbx.building.alternate_conformations.finalize_phil_str
}
""" % build_params_str

class rebuild_residue(object):
  """
  Callable wrapper class for rebuilding a single residue at a time.  This is
  not necessarily limited to modeling disorder, but it has been specifically
  designed to fit to a difference map in the presence of backbone and sidechain
  shifts.  Unlike some of the other tools, this method completely removes all
  sidechains within a sliding window before doing anything else.

  Only the target residue is returned; splitting of adjacent residues will be
  essential in most cases but is not handled here.
  """
  def __init__(self,
      target_map,
      pdb_hierarchy,
      xray_structure,
      geometry_restraints_manager,
      rotamer_eval,
      d_min):
    adopt_init_args(self, locals())
    from mmtbx.monomer_library import idealized_aa
    import mmtbx.monomer_library.server
    self.ideal_dict = idealized_aa.residue_dict()
    self.mon_lib_srv = mmtbx.monomer_library.server.server()

  def __call__(self,
      atom_group,
      log,
      window_size=2,
      backbone_sample_angle=10,
      anneal=False,
      annealing_temperature=1000,
      use_chi1_sampling=False):
    import iotbx.pdb.hierarchy
    from scitbx.array_family import flex
    assert (atom_group is not None)
    pdb_hierarchy = self.pdb_hierarchy.deep_copy()
    xray_structure = self.xray_structure.deep_copy_scatterers()
    geometry_restraints_manager = self.geometry_restraints_manager
    # FIXME this doesn't work - can't recover the atom_group afterwards!
    #hd_sel = xray_structure.hd_selection()
    #n_hydrogen = hd_sel.count(True)
    #if (n_hydrogen > 0):
    #  non_hd_sel = ~hd_sel
    #  pdb_hierarchy = pdb_hierarchy.select(non_hd_sel)
    #  xray_structure = xray_structure.select(non_hd_sel)
    #  geometry_restraints_manager = geometry_restraints_manager.select(
    #    non_hd_sel)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    isel = building.extract_iselection([atom_group])
    atom_group = pdb_atoms[isel[0]].parent()
    atom_group_start = atom_group.detached_copy()
    needs_rebuild = not building.is_stub_residue(atom_group)
    residue_group = atom_group.parent()
    assert (len(residue_group.atom_groups()) == 1)
    sel_residues = building.get_window_around_residue(
      residue=atom_group,
      window_size=window_size)
    # get rid of sidechains for surrounding residues only
    adjacent_residues = []
    for other_rg in sel_residues :
      if (other_rg != residue_group):
        adjacent_residues.append(other_rg)
    building.remove_sidechain_atoms(adjacent_residues)
    pdb_atoms = pdb_hierarchy.atoms()
    adjacent_trimmed_atom_names = pdb_atoms.extract_name()
    adjacent_trimmed_sel = pdb_atoms.extract_i_seq()
    xrs_adjacent_trimmed = xray_structure.select(adjacent_trimmed_sel)
    grm_adjacent_trimmed = geometry_restraints_manager.select(
      adjacent_trimmed_sel)
    pdb_atoms.reset_i_seq()
    # get rid of central sidechain and refine mainchain for entire window
    truncate = (not atom_group.resname in ["GLY","ALA"]) # XXX PRO?
    if (truncate):
      building.remove_sidechain_atoms([ atom_group ])
    pdb_atoms = pdb_hierarchy.atoms()
    all_mc_sel = pdb_atoms.extract_i_seq()
    xrs_mc = xrs_adjacent_trimmed.select(all_mc_sel)
    pdb_atoms.reset_i_seq()
    window_mc_sel = building.extract_iselection(sel_residues)
    selection = flex.bool(pdb_atoms.size(), False).set_selected(window_mc_sel,
      True)
    restraints_manager = grm_adjacent_trimmed.select(all_mc_sel)
    box = building.box_build_refine_base(
      xray_structure=xrs_mc,
      pdb_hierarchy=pdb_hierarchy,
      selection=selection,
      processed_pdb_file=None,
      target_map=self.target_map,
      geometry_restraints_manager=restraints_manager.geometry,
      d_min=self.d_min,
      out=null_out(),
      debug=True)
    box.restrain_atoms(
      selection=box.others_in_box,
      reference_sigma=0.1)
    box.real_space_refine(selection=box.selection_in_box)
    sites_new = box.update_original_coordinates()
    pdb_atoms.set_xyz(sites_new)
    # extend and replace existing residue.  this is done in such a way that
    # the original atom ordering for the central residue is preserved, which
    # allows us to use the pre-existing geometry restraints instead of
    # re-calculating them every time this function is called.
    target_atom_group = self.ideal_dict[atom_group.resname.lower()].\
      only_model().only_chain().only_residue_group().only_atom_group()
    new_atom_group_base = extend_sidechains.extend_residue(
      residue=atom_group,
      target_atom_group = target_atom_group,
      mon_lib_srv=self.mon_lib_srv)
    new_atom_group = iotbx.pdb.hierarchy.atom_group(resname=atom_group.resname)
    for atom in atom_group_start.atoms():
      for new_atom in new_atom_group_base.atoms():
        if (new_atom.name == atom.name):
          new_atom_group.append_atom(new_atom.detached_copy())
    n_atoms_new = len(new_atom_group.atoms())
    n_atoms_start = len(atom_group_start.atoms())
    if (n_atoms_new != n_atoms_start):
      raise RuntimeError(("Inconsistent atom counts for residue %s after "+
        "building (%d versus %d).") % (atom_group.id_str(), n_atoms_start,
        n_atoms_new))
    rg = atom_group.parent()
    rg.remove_atom_group(atom_group)
    rg.append_atom_group(new_atom_group)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    new_names = pdb_atoms.extract_name()
    assert new_names.all_eq(adjacent_trimmed_atom_names)
    # get new box around this residue
    residue_sel = building.extract_iselection([ new_atom_group ])
    selection = flex.bool(pdb_atoms.size(), False).set_selected(residue_sel,
      True)
    xrs_adjacent_trimmed.set_sites_cart(pdb_atoms.extract_xyz())
    box = building.box_build_refine_base(
      xray_structure=xrs_adjacent_trimmed,
      pdb_hierarchy=pdb_hierarchy,
      selection=selection,
      processed_pdb_file=None,
      target_map=self.target_map,
      geometry_restraints_manager=grm_adjacent_trimmed.geometry,
      d_min=self.d_min,
      out=null_out(),
      debug=True)
    # place sidechain using mmtbx.refinement.real_space.fit_residue
    if ((atom_group.resname in rotatable_sidechain_atoms) and
        (use_chi1_sampling)):
      fit_chi1_simple(
        residue=box.only_residue(),
        unit_cell=box.unit_cell_box,
        target_map=box.target_map_box,
        rotamer_eval=self.rotamer_eval)
      box.update_sites_from_pdb_atoms()
    else :
      box.fit_residue_in_box(backbone_sample_angle=backbone_sample_angle)
    if (anneal):
      box.anneal(start_temperature=annealing_temperature)
    #box.real_space_refine()
    sites_new = box.update_original_coordinates()
    pdb_hierarchy.atoms().set_xyz(sites_new)
    return building.atom_group_as_hierarchy(new_atom_group)

rotatable_sidechain_atoms = {
  "SER" : ["OG"],
  "CYS" : ["SG"],
}

def fit_chi1_simple(
    residue,
    unit_cell,
    target_map,
    rotamer_eval,
    sampling_angle=5):
  """
  For residues with only a single sidechain chi angle, we can instead sample
  the density at sidechain atoms directly instead of using the more opaque
  RSR infrastructure.  (Basically this is just the Ringer method applied to
  the Fo-Fc map.)  This is somewhat of a hack but should be more sensitive.
  """
  from scitbx.matrix import rotate_point_around_axis
  residue_atoms = residue.atoms()
  sites_start = residue_atoms.extract_xyz().deep_copy()
  sc_atoms = []
  c_alpha = c_beta = None
  assert (residue.resname in rotatable_sidechain_atoms)
  for atom in residue.atoms():
    name = atom.name.strip()
    if (name in rotatable_sidechain_atoms[residue.resname]):
      sc_atoms.append(atom)
    elif (name == "CA"):
      c_alpha = atom
    elif (name == "CB"):
      c_beta = atom
  if (len(sc_atoms) == 0) or (None in [c_alpha, c_beta]):
    return False
  angle = 0
  map_level_best = 3.0 * len(sc_atoms)
  sites_best = sites_start
  is_acceptable_cys_sg = True
  while (angle < 360):
    map_level_sum = 0
    for atom in sc_atoms :
      atom.xyz = rotate_point_around_axis(
        axis_point_1=c_alpha.xyz,
        axis_point_2=c_beta.xyz,
        point=atom.xyz,
        angle=sampling_angle,
        deg=True)
      site_frac = unit_cell.fractionalize(site_cart=atom.xyz)
      map_level = target_map.eight_point_interpolation(site_frac)
      map_level_sum += map_level
      if (residue.resname == "CYS") and (atom.name.strip() == "SG"):
        if (map_level < 5.0):
          is_acceptable_cys_sg = False
    angle += sampling_angle
    if (not is_acceptable_cys_sg) : # always True if resname != CYS
      continue
    rotamer = rotamer_eval.evaluate_residue(residue)
    #print angle, map_level_sum, rotamer
    if (map_level_sum > map_level_best) and (rotamer != "OUTLIER"):
      #print "setting sites_best"
      sites_best = residue_atoms.extract_xyz().deep_copy()
      map_level_best = map_level_sum
  #print "RMSD:", sites_best.rms_difference(sites_start)
  residue_atoms.set_xyz(sites_best)

class residue_trial(slots_getstate_setstate):
  __slots__ = [ "new_hierarchy", "sc_n_atoms", "sc_two_fofc_mean",
                "sc_fofc_mean", "two_fofc_values", "fofc_values",
                "stats", "occupancy", "rotamer", ]
  def __init__(self, residue, new_hierarchy, occupancy, rotamer_eval,
                      fmodel, two_fofc_map, fofc_map):
    self.new_hierarchy = new_hierarchy
    self.occupancy = occupancy
    self.two_fofc_values = []
    self.fofc_values = []
    self.sc_n_atoms = 0
    self.sc_fofc_mean = self.sc_two_fofc_mean = None
    unit_cell = fmodel.xray_structure.unit_cell()
    sc_fofc_sum = sc_two_fofc_sum = 0
    for atom in new_hierarchy.atoms():
      assert (not atom.element.strip() in ["H","D"])
      name = atom.name.strip()
      site_frac = unit_cell.fractionalize(site_cart=atom.xyz)
      two_fofc_value = two_fofc_map.eight_point_interpolation(site_frac)
      fofc_value = fofc_map.eight_point_interpolation(site_frac)
      self.two_fofc_values.append(two_fofc_value)
      self.fofc_values.append(fofc_value)
      if (not name in ["N","C","O","CA", "CB"]):
        self.sc_n_atoms += 1
        sc_fofc_sum += fofc_value
        sc_two_fofc_sum += two_fofc_value
    if (self.sc_n_atoms > 0):
      self.sc_fofc_mean = sc_fofc_sum / self.sc_n_atoms
      self.sc_two_fofc_mean = sc_two_fofc_sum / self.sc_n_atoms
    new_conf = self.as_atom_group()
    self.stats = alt_confs.coord_stats_for_atom_groups(residue, new_conf)
    self.rotamer = None
    if (not residue.resname in ["GLY","ALA","PRO"]):
      self.rotamer = rotamer_eval.evaluate_residue(new_conf)

  def as_atom_group(self):
    return self.new_hierarchy.only_model().only_chain().only_residue_group().\
      only_atom_group()

  def rescore(self, params, log=None, prefix=""):
    # FIXME this needs to be rationalized and made consistent with the
    # rescoring against the original maps
    if (log is None) : log = null_out()
    reject = (self.rotamer == "OUTLIER")
    bad_mc_two_fofc_msg = None
    for i_seq, atom in enumerate(self.new_hierarchy.atoms()):
      name = atom.name.strip()
      if (name in ["N","C","O","CA", "CB"]):
        if (self.two_fofc_values[i_seq] < params.omit_two_fofc_min_mc):
          bad_mc_two_fofc_msg = "poor backbone: 2Fo-Fc(%s)=%.2f" % (name,
                self.two_fofc_values[i_seq])
          reject = True
          break
    if (self.sc_fofc_mean is not None):
      if (self.sc_fofc_mean < params.omit_fofc_min_sc_mean):
        reject = True
      elif ((self.sc_two_fofc_mean < params.omit_two_fofc_min_sc_mean) and
            (self.sc_fofc_mean < params.omit_fofc_min_sc_mean + 1)):
        reject = True
    flag = ""
    if (reject):
      flag = " !!!"
    print(prefix+\
      "occupancy=%.2f rotamer=%s 2Fo-Fc(sc)=%s  Fo-Fc(sc)=%s%s" % \
      (self.occupancy, self.rotamer,
       format_value("%5f", self.sc_two_fofc_mean),
       format_value("%5f", self.sc_fofc_mean),
       flag), file=log)
    if (bad_mc_two_fofc_msg is not None):
      print(prefix+"  %s" % bad_mc_two_fofc_msg, file=log)
    return (not reject)

def find_alternate_residue(residue,
    pdb_hierarchy,
    fmodel,
    restraints_manager,
    params,
    verbose=False,
    debug=None,
    log=None):
  if (log is None):
    log = null_out()
  t1 = time.time()
  from scitbx.array_family import flex
  selection = flex.size_t()
  window = building.get_window_around_residue(residue,
    window_size=params.window_size)
  for pdb_object in window :
    selection.extend(pdb_object.atoms().extract_i_seq())
  assert (len(selection) > 0) and (not selection.all_eq(0))
  occupancies = []
  if (params.expected_occupancy is not None):
    assert (0.0 <= params.expected_occupancy <= 1.0)
    occupancies = [ params.expected_occupancy ]
  else :
    occupancies = [ 0.2, 0.3, 0.4, 0.5 ]
  trials = []
  sites_start_1d = pdb_hierarchy.atoms().extract_xyz().as_double()
  from mmtbx.rotamer import rotamer_eval
  rotamer_manager = rotamer_eval.RotamerEval(data_version="8000")
  id_str = residue.id_str()
  delete_selection = None
  if (params.omit_waters):
    delete_selection = building.get_nearby_water_selection(
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=fmodel.xray_structure,
      selection=selection)
  for occupancy in occupancies :
    prefix = "%s_%.2f" % (id_str.replace(" ", "_"), occupancy)
    map_file_name = None
    if (debug > 1):
      map_file_name = prefix + ".mtz"
    two_fofc_map, fofc_map = alt_confs.get_partial_omit_map(
      fmodel=fmodel.deep_copy(),
      selection=selection,
      selection_delete=delete_selection,
      negate_surrounding=True,
      map_file_name=map_file_name,
      partial_occupancy=1.0 - occupancy)
    rebuild = rebuild_residue(
      target_map=fofc_map,
      pdb_hierarchy=pdb_hierarchy,
      xray_structure=fmodel.xray_structure,
      geometry_restraints_manager=restraints_manager,
      rotamer_eval=rotamer_manager,
      d_min=fmodel.f_obs().d_min())
    new_hierarchy = rebuild(atom_group=residue,
      window_size=params.window_size,
      backbone_sample_angle=params.backbone_sample_angle,
      anneal=params.anneal,
      annealing_temperature=params.annealing_temperature,
      use_chi1_sampling=params.simple_chi1_sampling,
      log=log)
    trial = residue_trial(
      residue=residue,
      new_hierarchy=new_hierarchy,
      occupancy=occupancy,
      rotamer_eval=rotamer_manager,
      fmodel=fmodel,
      two_fofc_map=two_fofc_map,
      fofc_map=fofc_map)
    trials.append(trial)
    if (debug > 1):
      open("%s.pdb" % prefix, "w").write(trial.new_hierarchy.as_pdb_string())
  sites_end_1d = pdb_hierarchy.atoms().extract_xyz().as_double()
  assert sites_start_1d.all_eq(sites_end_1d)
  t2 = time.time()
  if (debug > 1):
    print("  %d build trials (%s): %.3fs" % (len(occupancies),
      residue.id_str(), t2 - t1), file=log)
  return trials

class find_all_alternates(object):
  """
  Wrapper for parallelizing calls to find_alternate_residue.
  """
  def __init__(self,
      residues,
      pdb_hierarchy,
      fmodel,
      restraints_manager,
      params,
      nproc=Auto,
      verbose=False,
      debug=None,
      log=sys.stdout):
    adopt_init_args(self, locals())
    nproc = easy_mp.get_processes(nproc)
    print("", file=log)
    if (nproc == 1):
      print("  running all residues serially", file=log)
      self.results = []
      for i_res in range(len(residues)):
        self.results.append(self.__call__(i_res, log=log))
    else :
      print("  will use %d processes" % nproc, file=log)
      self.results = easy_mp.pool_map(
        fixed_func=self,
        iterable=range(len(residues)),
        processes=nproc)

  def __call__(self, i_res, log=None):
    return find_alternate_residue(
      residue=self.residues[i_res],
      pdb_hierarchy=self.pdb_hierarchy,
      fmodel=self.fmodel,
      restraints_manager=self.restraints_manager,
      params=self.params,
      verbose=self.verbose,
      debug=self.debug,
      log=log)

def pick_best_alternate(
    trials,
    params,
    rotamer,
    log=None):
  if (log is None):
    log = null_out()
  filtered = []
  for trial in trials :
    accept = trial.rescore(params=params.map_thresholds,
      log=log,
      prefix="    ")
    if (accept) : filtered.append(trial)
  trials = filtered
  if (len(trials) == 0):
    return None
  elif (len(trials) == 1):
    return trials[0]
  else :
    rmsd_max = 0
    best_trial = None
    for trial in trials :
      if (trial.stats.rmsd > rmsd_max):
        rmsd_max = trial.stats.rmsd
        best_trial = trial
    return best_trial

def process_results(
    pdb_hierarchy,
    fmodel,
    residues_in,
    building_trials,
    params,
    verbose=False,
    log=sys.stdout):
  assert (len(residues_in) == len(building_trials))
  from mmtbx.rotamer import rotamer_eval
  n_alternates = 0
  unit_cell = fmodel.xray_structure.unit_cell()
  two_fofc_map, fofc_map = building.get_difference_maps(fmodel)
  rot_eval = rotamer_eval.RotamerEval(data_version="8000")
  for main_conf, trials in zip(residues_in, building_trials):
    if (trials is None):
      print("WARNING: error building %s" % main_conf.id_str(), file=log)
      continue
    if (len(trials) == 0):
      continue
    res_log = StringIO()
    print("  %s:" % main_conf.id_str(), file=res_log)
    main_rotamer = alt_rotamer = None
    if (not main_conf.resname in ["GLY","PRO","ALA"]):
      main_rotamer = rot_eval.evaluate_residue(main_conf)
      assert (main_rotamer != "OUTLIER")
    best_trial = pick_best_alternate(
      trials=trials,
      params=params,
      rotamer=main_rotamer,
      log=res_log)
    if (best_trial is None):
      continue
    new_conf = best_trial.as_atom_group()
    changed_rotamer = (best_trial.rotamer != main_rotamer)
    skip = False
    flag = ""
    stats = best_trial.stats
    # FIXME this needs to be made more consistent with the filtering criteria
    # in disorder/__init__.py
    if ((stats.rmsd < params.rmsd_min) and
        (stats.max_dev < params.rmsd_min) and
        (not changed_rotamer)):
      skip = True
    print("    selected conformer (occ=%.2f):" % \
      best_trial.occupancy, file=res_log)
    res_log2 = StringIO()
    density_quality = building.residue_density_quality(
      atom_group=new_conf,
      unit_cell=unit_cell,
      two_fofc_map=two_fofc_map,
      fofc_map=fofc_map)
    if (main_conf.resname in ["CYS", "MET", "MSE"]):
      # XXX if we have a heavier atom (S or SE) in the residue, some additional
      # sanity checks insure that it has slightly stronger density than we
      # require for lighter elements.  multipliers are just guesses, taking
      # into account rad damage and partial SeMet incorporation.
      heavy_atom = {"CYS":"SG", "MET":"SD", "MSE":"SE"}[main_conf.resname]
      mult = {"CYS":1.6, "MET":1.6, "MSE":2.0}[main_conf.resname]
      map_levels = density_quality.density_at_atom(heavy_atom)
      if (map_levels is None) : # this probably shouldn't even happen
        skip = True
      if ((map_levels.fofc < params.map_thresholds.fofc_min*mult) or
          (map_levels.two_fofc < params.map_thresholds.two_fofc_min*mult)):
        skip = True
    n_atoms_outside_density = density_quality.show_atoms_outside_density(
      two_fofc_cutoff=params.map_thresholds.two_fofc_min,
      fofc_cutoff=params.map_thresholds.fofc_min,
      out=res_log2,
      prefix="      ")
    fofc_max = density_quality.max_fofc_value()
    if (n_atoms_outside_density > 0):
      skip = True
    if (not skip) and (verbose):
      flag = " ***"
    print("      RMSD=%5.3f  max. change=%.2f  max(Fo-Fc)=%.1f%s" \
      % (stats.rmsd, stats.max_dev, fofc_max, flag), file=res_log)
    if (changed_rotamer):
      print("      starting rotamer=%s  new rotamer=%s" % \
        (main_rotamer, best_trial.rotamer), file=res_log)
    if (n_atoms_outside_density != 0):
      print("      atoms outside density:", file=res_log)
      print(res_log2.getvalue(), file=res_log)
    else :
      print("", file=res_log)
    if (not skip) or (verbose):
      log.write(res_log.getvalue())
    if (skip) : continue
    residue_group = main_conf.parent()
    main_conf.altloc = 'A'
    new_occ = 0.5
    if (params.expected_occupancy is not None):
      new_occ = max(0.2, min(0.8, params.expected_occupancy))
    for atom in main_conf.atoms():
      atom.occ = 1.0 - new_occ
      atom.segid = alt_confs.SEGID_MAIN
    new_conf = new_conf.detached_copy()
    new_conf.altloc = 'B'
    for atom in new_conf.atoms():
      atom.segid = alt_confs.SEGID_NEW_REBUILT
      atom.occ = new_occ
    residue_group.append_atom_group(new_conf)
    n_alternates += 1
  return n_alternates

class rsr_fragments_parallel(alt_confs.rsr_fragments_base):
  def __call__(self, fragment):
    from scitbx.array_family import flex
    frag_selection = flex.bool(self.pdb_atoms.size(), fragment)
    sites_cart_orig = self.fmodel.xray_structure.sites_cart()
    sites_start = sites_cart_orig.select(fragment)
    atom_start = self.pdb_atoms[fragment[0]].fetch_labels()
    atom_end = self.pdb_atoms[fragment[-1]].fetch_labels()
    label = "%s%s-%s" % (atom_start.chain_id, atom_start.resid(),
      atom_end.resid())
    # XXX should the partial occupancy be determined automatically?
    two_fofc_map, fofc_map = alt_confs.get_partial_omit_map(
      fmodel=self.fmodel.deep_copy(),
      selection=(frag_selection & self.sele_main_conf),
      selection_delete=(frag_selection & ~(self.sele_main_conf)),
      partial_occupancy=0.6)
    target_map = two_fofc_map
    if (self.rsr_fofc_map_target):
      target_map = fofc_map
    box = self.make_box(selection=frag_selection, target_map=target_map)
    box.restrain_atoms(
      selection=alt_confs.SELECTION_OLD,
      reference_sigma=0.02,
      reset_first=True)
    box.restrain_atoms(
      selection=alt_confs.SELECTION_NEW_REBUILT,
      reference_sigma=0.05,
      reset_first=False)
    # first fix the geometry of adjacent residues
    box.real_space_refine(alt_confs.SELECTION_NEW_SPLIT)
    # now the entire B conformer
    box.real_space_refine(alt_confs.SELECTION_NEW)
    sites_cart_refined = box.update_original_coordinates()
    sites_new = sites_cart_refined.select(fragment)
    self.fmodel.xray_structure.set_sites_cart(sites_cart_orig)
    self.pdb_atoms.set_xyz(sites_cart_orig)
    return alt_confs.refined_fragment(
      label=label,
      selection=frag_selection,
      sites_cart=sites_cart_refined,
      rmsd=sites_new.rms_difference(sites_start))

def real_space_refine(
    pdb_hierarchy,
    fmodel,
    cif_objects,
    params,
    out,
    nproc=None,
    max_cycles=100, # arbitrarily large
    remediate=False):
  from scitbx.array_family import flex
  i_cycle = 0
  while (i_cycle < max_cycles):
    print("  Cycle %d:" % (i_cycle+1), file=out)
    # this keeps track of which residues were split in the previous cycle -
    # we only refine segments that have had residues added
    rebuilt_flags = pdb_hierarchy.atoms().extract_tmp_as_size_t()
    processed_pdb_file = building.reprocess_pdb(
      pdb_hierarchy=pdb_hierarchy,
      cif_objects=cif_objects,
      crystal_symmetry=fmodel.xray_structure,
      out=null_out())
    # get the 2mFo-DFc map without new alternates!
    #two_fofc_map = fmodel.two_fofc_map(exclude_free_r_reflections=True)
    pdb_hierarchy = processed_pdb_file.all_chain_proxies.pdb_hierarchy
    pdb_atoms = pdb_hierarchy.atoms()
    xray_structure = processed_pdb_file.xray_structure()
    geometry_restraints_manager = \
      processed_pdb_file.geometry_restraints_manager(show_energies=False)
    fmodel.update_xray_structure(xray_structure)
    sele_cache = pdb_hierarchy.atom_selection_cache()
    # FIXME very inefficient when looping!
    # this will include both the newly built residues and the original atoms,
    # including residues split to allow for backbone flexibility.
    sele_split = sele_cache.selection(alt_confs.SELECTION_MODIFIED)
    sele_main_conf = sele_cache.selection(alt_confs.SELECTION_OLD)
    assert (len(sele_split) > 0)
    k = 0
    fragments = []
    while (k < len(sele_split)):
      if (sele_split[k]):
        current_fragment = flex.size_t()
        while (sele_split[k]):
          current_fragment.append(k)
          k += 1
        atom_start = pdb_atoms[current_fragment[0]].fetch_labels()
        atom_end = pdb_atoms[current_fragment[-1]].fetch_labels()
        frag_selection = flex.bool(sele_split.size(), current_fragment)
        if (i_cycle > 0):
          flags = rebuilt_flags.select(frag_selection)
          if flags.all_eq(0):
            continue
        fragments.append(current_fragment)
      else :
        k += 1
    if (len(fragments) == 0):
      pass
    refine_fragments = rsr_fragments_parallel(
      pdb_hierarchy=pdb_hierarchy,
      fmodel=fmodel,
      processed_pdb_file=processed_pdb_file,
      sele_main_conf=sele_main_conf,
      rsr_fofc_map_target=(i_cycle==0 and params.cleanup.rsr_fofc_map_target))
    refined = easy_mp.pool_map(
      fixed_func=refine_fragments,
      iterable=fragments,
      processes=nproc)
    sites_refined = pdb_atoms.extract_xyz()
    for result in refined :
      assert (result is not None)
      result.show(out=out, prefix="    ")
      sites_refined.set_selected(result.selection, result.sites_cart)
    pdb_atoms.set_xyz(sites_refined)
    xray_structure.set_sites_cart(sites_refined)
    fmodel.update_xray_structure(xray_structure)
    if (not remediate) or (max_cycles == 1):
      break
    else :
      for atom in pdb_hierarchy.atoms():
        if (atom.segid == alt_confs.SEGID_NEW_SPLIT):
          atom.segid = alt_confs.SEGID_NEW_REBUILT
      print("    checking for conformational strain...", file=out)
      n_split = alt_confs.spread_alternates(
        pdb_hierarchy=pdb_hierarchy,
        new_occupancy=params.residue_fitting.expected_occupancy,
        split_all_adjacent=False,
        selection=alt_confs.SELECTION_NEW_REBUILT)
      if (n_split > 0):
        print("    split another %d residue(s) - will re-run RSR" % \
          n_split, file=out)
      else :
        break
    i_cycle += 1
  xray_structure = pdb_hierarchy.extract_xray_structure(
    crystal_symmetry=fmodel.xray_structure)
  fmodel.update_xray_structure(xray_structure,
    update_f_mask=True,
    update_f_calc=True)
  #fmodel.info().show_targets(out=out, text="After real-space refinement")
  t2 = time.time()
  return pdb_hierarchy

def build_cycle(pdb_hierarchy,
    fmodel,
    geometry_restraints_manager,
    params,
    selection=None,
    cif_objects=(),
    nproc=Auto,
    out=sys.stdout,
    verbose=False,
    debug=None,
    i_cycle=0):
  from mmtbx import restraints
  from scitbx.array_family import flex
  t_start = time.time()
  hd_sel = fmodel.xray_structure.hd_selection()
  n_hydrogen = hd_sel.count(True)
  if (n_hydrogen > 0) and (True) : #params.building.delete_hydrogens):
    print("WARNING: %d hydrogen atoms will be removed!" % n_hydrogen, file=out)
    non_hd_sel = ~hd_sel
    # XXX it's better to do this in-place for the hierarchy, because calling
    # pdb_hierarchy.select(non_hd_sel) will not remove parent-child
    # relationships involving hydrogens, which causes problems when running
    # the MolProbity validation.
    pdb_hierarchy.remove_hd(reset_i_seq=True)
    xray_structure = fmodel.xray_structure.select(non_hd_sel)
    assert (pdb_hierarchy.atoms_size() == xray_structure.scatterers().size())
    fmodel.update_xray_structure(xray_structure)
    geometry_restraints_manager = geometry_restraints_manager.select(non_hd_sel)
  pdb_atoms = pdb_hierarchy.atoms()
  segids = pdb_atoms.extract_segid().strip()
  if (not segids.all_eq("")):
    print("WARNING: resetting segids to blank", file=out)
    for i_seq, atom in enumerate(pdb_atoms):
      atom.segid = ""
      sc = fmodel.xray_structure.scatterers()[i_seq]
      sc.label = atom.id_str()
  if isinstance(selection, str):
    sele_cache = pdb_hierarchy.atom_selection_cache()
    selection = sele_cache.selection(selection)
  make_header("Build cycle %d" % (i_cycle+1), out=out)
  fmodel.info().show_rfactors_targets_scales_overall(out=out)
  if (debug > 0):
    from mmtbx.maps.utils import get_maps_from_fmodel
    from iotbx.map_tools import write_map_coeffs
    two_fofc, fofc = get_maps_from_fmodel(fmodel,
      exclude_free_r_reflections=True)
    write_map_coeffs(
      fwt_coeffs=two_fofc,
      delfwt_coeffs=fofc,
      file_name="cycle_%d_start.mtz" % (i_cycle+1))
  candidate_residues = alt_confs.filter_before_build(
    pdb_hierarchy=pdb_hierarchy,
    fmodel=fmodel,
    geometry_restraints_manager=geometry_restraints_manager,
    selection=selection,
    params=params.prefilter,
    verbose=verbose,
    log=out)
  t1 = time.time()
  print("filtering: %.3fs" % (t1-t_start), file=out)
  restraints_manager = restraints.manager(
    geometry=geometry_restraints_manager,
    normalization=True)
  make_sub_header("Finding alternate conformations", out=out)
  building_trials = find_all_alternates(
    residues=candidate_residues,
    pdb_hierarchy=pdb_hierarchy,
    restraints_manager=restraints_manager,
    fmodel=fmodel,
    params=params.residue_fitting,
    nproc=params.nproc,
    verbose=verbose,
    debug=debug,
    log=out).results
  t2 = time.time()
  print("  building: %.3fs" % (t2-t1), file=out)
  make_sub_header("Scoring and assembling alternates", out=out)
  n_alternates = process_results(
    pdb_hierarchy=pdb_hierarchy,
    fmodel=fmodel,
    residues_in=candidate_residues,
    building_trials=building_trials,
    params=params.residue_fitting,
    verbose=verbose,
    log=out)
  if (n_alternates > 0):
    print("", file=out)
    print("  %d disordered residues built" % n_alternates, file=out)
    n_split = alt_confs.spread_alternates(pdb_hierarchy,
      new_occupancy=params.residue_fitting.expected_occupancy,
      split_all_adjacent=True,
      log=out)
    assert (n_split > 0)
    print("  %d adjacent residues split" % n_split, file=out)
  else :
    print("No alternates built this round.", file=out)
  t3 = time.time()
  print("  assembly: %.3fs" % (t3-t2), file=out)
  if (not params.cleanup.rsr_after_build):
    if (n_alternates > 0):
      print("Skipping final RSR step (rsr_after_build=False).", file=out)
    else :
      print("No refinement needs to be performed.", file=out)
  else :
    make_sub_header("Real-space refinement", out=out)
    print("", file=out)
    pdb_hierarchy = real_space_refine(
      pdb_hierarchy=pdb_hierarchy,
      fmodel=fmodel,
      cif_objects=cif_objects,
      params=params,
      nproc=params.nproc,
      remediate=True,
      out=out)
    t4 = time.time()
    print("", file=out)
    print("RSR: %.3fs" % (t4-t3), file=out)
  fmodel.info().show_targets(out=out, text="Rebuilt model")
  t_end = time.time()
  alt_confs.finalize_model(
    pdb_hierarchy=pdb_hierarchy,
    xray_structure=pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=fmodel.xray_structure),
    set_b_iso=params.cleanup.set_b_iso,
    convert_to_isotropic=params.cleanup.convert_to_isotropic,
    selection="altloc A or altloc B")
  t_end = time.time()
  print("Total runtime for cycle: %.3fs" % (t_end-t_start), file=out)
  return pdb_hierarchy, n_alternates
