
# TODO avoid re-processing structure each time
# TODO tests

"""
Prototype for building alternate conformations into difference density.
The actual method is a variant of one that Pavel suggested, in combination
with the procedure in mmtbx.building.extend_sidechains: first, the
backbone atoms for a residue and its neighbors are refined into the
target map using minimization and/or annealing, then the sidechain is
replaced (using an idealized copy) and its placement optimized by a grid
search that also allows for backbone flexibility.
"""

from __future__ import division
from mmtbx.building import extend_sidechains
from mmtbx.building import disorder
from mmtbx import building
from libtbx import adopt_init_args, Auto
from libtbx.utils import null_out
from libtbx import easy_mp
import sys

class rebuild_residue (object) :
  """
  Callable wrapper class for rebuilding a single residue at a time.  This is
  not necessarily limited to modeling disorder, but it has been specifically
  designed to fit to a difference map in the presence of backbone and sidechain
  shifts.  Unlike some of the other tools, this method completely removes all
  sidechains within a sliding window before doing anything else.

  Only the target residue is returned; splitting of adjacent residues will be
  essential in most cases but is not handled here.
  """
  def __init__ (self,
      target_map,
      pdb_hierarchy,
      xray_structure,
      geometry_restraints_manager,
      d_min) :
    adopt_init_args(self, locals())
    from mmtbx.monomer_library import idealized_aa
    from mmtbx.rotamer import rotamer_eval
    import mmtbx.monomer_library.server
    self.ideal_dict = idealized_aa.residue_dict()
    self.mon_lib_srv = mmtbx.monomer_library.server.server()
    self.rotamer_manager = rotamer_eval.RotamerEval()

  def __call__ (self, atom_group, log, window_size=2, anneal=False) :
    from scitbx.array_family import flex
    assert (atom_group is not None)
    pdb_hierarchy = self.pdb_hierarchy.deep_copy()
    xray_structure = self.xray_structure.deep_copy_scatterers()
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    isel = building.extract_iselection([atom_group])
    atom_group = pdb_atoms[isel[0]].parent()
    residue_group = atom_group.parent()
    assert (len(residue_group.atom_groups()) == 1)
    sel_residues = building.get_window_around_residue(
      residue=atom_group,
      window_size=window_size)
    print >> log, "  %d residues extracted" % len(sel_residues)
    print >> log, "  removing sidechain atoms..."
    building.remove_sidechain_atoms(sel_residues)
    pdb_atoms = pdb_hierarchy.atoms()
    all_mc_sel = pdb_atoms.extract_i_seq()
    xrs_mc = xray_structure.select(all_mc_sel)
    pdb_atoms.reset_i_seq()
    window_mc_sel = building.extract_iselection(sel_residues)
    selection = flex.bool(pdb_atoms.size(), False).set_selected(window_mc_sel,
      True)
    restraints_manager = self.geometry_restraints_manager.select(all_mc_sel)
    box = building.box_build_refine_base(
      xray_structure=xrs_mc,
      pdb_hierarchy=pdb_hierarchy,
      selection=selection,
      processed_pdb_file=None,
      target_map=self.target_map,
      geometry_restraints_manager=restraints_manager.geometry,
      d_min=self.d_min,
      out=log,
      debug=True)
    box.restrain_atoms(
      selection=box.others_in_box,
      reference_sigma=0.1)
    box.real_space_refine(selection=box.selection_in_box)
    if (anneal) : # TODO
      raise NotImplementedError()
    sites_new = box.update_original_coordinates()
    pdb_atoms.set_xyz(sites_new)
    # extend and replace existing residue
    new_atom_group = extend_sidechains.extend_residue(
      residue=atom_group,
      ideal_dict=self.ideal_dict,
      hydrogens=False,
      mon_lib_srv=self.mon_lib_srv,
      match_conformation=True)
    rg = atom_group.parent()
    rg.remove_atom_group(atom_group)
    rg.append_atom_group(new_atom_group)
    pdb_atoms = pdb_hierarchy.atoms()
    pdb_atoms.reset_i_seq()
    # get new box around this residue
    residue_sel = building.extract_iselection([ new_atom_group ])
    selection = flex.bool(pdb_atoms.size(), False).set_selected(residue_sel,
      True)
    xray_structure = pdb_hierarchy.extract_xray_structure(
      crystal_symmetry=self.xray_structure)
    # FIXME this is horrendously inefficient
    box = building.box_build_refine_base(
      xray_structure=xray_structure,
      pdb_hierarchy=pdb_hierarchy,
      selection=selection,
      processed_pdb_file=None,
      target_map=self.target_map,
      d_min=self.d_min,
      out=null_out(),
      debug=True)
    box.fit_residue_in_box(backbone_sample_angle=5)
    sites_new = box.update_original_coordinates()
    pdb_hierarchy.atoms().set_xyz(sites_new)
    return building.atom_group_as_hierarchy(new_atom_group)

params_str = """
  expected_occupancy = 0.4
    .type = float
  window_size = 2
    .type = int
  rescore = True
    .type = bool
  map_thresholds {
    two_fofc_min_sc_mean = 0.8
      .type = float
    two_fofc_min_mc = 1.0
      .type = float
    fofc_min_sc_mean = 2.5
      .type = float
    starting_fofc_min_sc_single = 3.0
      .type = float
  }
"""

def find_alternate_residue (residue,
    pdb_hierarchy,
    fmodel,
    restraints_manager,
    params,
    log=None) :
  if (log is None) :
    log = null_out()
  from scitbx.array_family import flex
  selection = flex.size_t()
  window = building.get_window_around_residue(residue,
    window_size=params.window_size)
  for pdb_object in window :
    selection.extend(pdb_object.atoms().extract_i_seq())
  assert (len(selection) > 0) and (not selection.all_eq(0))
  two_fofc_map, fofc_map = disorder.get_partial_omit_map(
    fmodel=fmodel,
    selection=selection,
    selection_delete=None,#nearby_water_selection,
    negate_surrounding=True,
    partial_occupancy=1.0 - params.expected_occupancy)
  rebuild = rebuild_residue(
    target_map=fofc_map,
    pdb_hierarchy=pdb_hierarchy,
    xray_structure=fmodel.xray_structure,
    geometry_restraints_manager=restraints_manager,
    d_min=fmodel.f_obs().d_min())
  new_hierarchy = rebuild(atom_group=residue,
    window_size=params.window_size,
    log=log)
  if (params.rescore) :
    unit_cell = fmodel.xray_structure.unit_cell()
    sc_total = sc_fofc_sum = sc_two_fofc_sum = 0
    for atom in new_hierarchy.atoms() :
      name = atom.name.strip()
      site_frac = unit_cell.fractionalize(site_cart=atom.xyz)
      two_fofc_value = two_fofc_map.eight_point_interpolation(site_frac)
      fofc_value = fofc_map.eight_point_interpolation(site_frac)
      if (name in ["N","C","O","CA", "CB"]) :
        if (two_fofc_value < params.map_thresholds.two_fofc_min_mc) :
          return None
      else :
        sc_total += 1
        sc_fofc_sum += fofc_value
        sc_two_fofc_sum += two_fofc_value
    if (sc_total > 0) :
      sc_fofc_mean = sc_fofc_sum / sc_total
      sc_two_fofc_mean = sc_two_fofc_sum / sc_total
      if (sc_fofc_mean < params.map_thresholds.fofc_min_sc_mean) :
        return None
      if (sc_two_fofc_mean < params.map_thresholds.two_fofc_min_sc_mean) :
        return None
  return new_hierarchy

class find_all_alternates (object) :
  """
  Wrapper for parallelizing calls to find_alternate_residue.
  """
  def __init__ (self,
      residues,
      pdb_hierarchy,
      fmodel,
      restraints_manager,
      params,
      nproc=Auto,
      log=sys.stdout) :
    adopt_init_args(self, locals())
    nproc = easy_mp.get_processes(nproc)
    print >> log, "Will use %d processes" % nproc
    if (nproc == 1) :
      self.results = []
      for i_res in range(len(residues)) :
        self.results.append(self.__call__(i_res, log=log))
    else :
      self.results = easy_mp.pool_map(
        fixed_func=self,
        iterable=range(len(residues)),
        processes=nproc)

  def __call__ (self, i_res, log=None) :
    return find_alternate_residue(
      residue=self.residues[i_res],
      pdb_hierarchy=self.pdb_hierarchy,
      fmodel=self.fmodel,
      restraints_manager=self.restraints_manager,
      params=self.params,
      log=log)

def process_results (
    pdb_hierarchy,
    fmodel,
    residues_in,
    new_residues,
    params,
    log=sys.stdout) :
  assert (len(residues_in) == len(new_residues))
  residues_out = []
  for h in new_residues :
    if (h is not None) :
      residues_out.append(
        h.only_model().only_chain().only_residue_group().only_atom_group())
    else :
      residues_out.append(None)
  unit_cell = fmodel.xray_structure.unit_cell()
  two_fofc_map = fofc_map = None
  if (params.rescore) :
    two_fofc_map, fofc_map = building.get_difference_maps(fmodel)
  print >> log, ""
  print >> log, "Assembling disordered residues..."
  for main_conf, new_conf in zip(residues_in, residues_out) :
    if (new_conf is None) : continue
    stats = disorder.coord_stats_for_atom_groups(main_conf, new_conf)
    if (stats.rmsd > 0.5) :
      if (params.rescore) :
        fofc_max = 0
        for atom in new_conf.atoms() :
          name = atom.name.strip()
          site_frac = unit_cell.fractionalize(site_cart=atom.xyz)
          fofc_value = fofc_map.eight_point_interpolation(site_frac)
          if ((not name in ["N","C","O","CA", "CB"]) or
              ((name == "O") and (new_conf.resname in ["GLY", "ALA"]))) :
            if (fofc_value > fofc_max) :
              fofc_max = fofc_value
        if (fofc_max < params.map_thresholds.starting_fofc_min_sc_single) :
          continue
        print >> log, "  %s: RMSD=%5.3f  max. change=%.2f  max(Fo-Fc)=%.1f" % \
          (main_conf.id_str(), stats.rmsd, stats.max_dev, fofc_max)
      else :
        print >> log, "  %s: RMSD=%5.3f  max. change=%.2f" % \
          (main_conf.id_str(), stats.rmsd, stats.max_dev)
      residue_group = main_conf.parent()
      main_conf.altloc = 'A'
      for atom in main_conf.atoms() :
        atom.segid = "XXXX"
        atom.occ = 1.0 - params.expected_occupancy
      alt_conf = new_conf.detached_copy()
      alt_conf.altloc = 'B'
      for atom in alt_conf.atoms() :
        atom.segid = "XXXX"
        atom.occ = params.expected_occupancy
      residue_group.append_atom_group(alt_conf)
  spread_alternates(pdb_hierarchy, params, log=log)

def spread_alternates (pdb_hierarchy, params, log) :
  print >> log, ""
  print >> log, "Splitting adjacent residues..."
  def split_residue (residue_group) :
    print >> log, "  %s" % residue_group.id_str()
    main_conf = residue_group.only_atom_group()
    main_conf.altloc = 'A'
    for atom in main_conf.atoms() :
      atom.occ = 1.0 - params.expected_occupancy
      atom.segid = 'ZZZZ'
    alt_conf = main_conf.detached_copy()
    alt_conf.altloc = 'B'
    for atom in alt_conf.atoms() :
      atom.occ = params.expected_occupancy
      atom.segid = 'ZZZZ'
    residue_group.append_atom_group(alt_conf)
  for chain in pdb_hierarchy.only_model().chains() :
    residue_groups = chain.residue_groups()
    for i_res, residue_group in enumerate(residue_groups) :
      atom_groups = residue_group.atom_groups()
      if (len(atom_groups) > 1) :
        segid = residue_group.atoms()[0].segid
        if (segid != 'XXXX') :
          continue
        if (i_res > 0) :
          prev_group = residue_groups[i_res - 1]
          if (len(prev_group.atom_groups()) == 1) :
            split_residue(prev_group)
        if (i_res < len(residue_groups) - 1) :
          next_group = residue_groups[i_res + 1]
          if (len(next_group.atom_groups()) == 1) :
            split_residue(next_group)
