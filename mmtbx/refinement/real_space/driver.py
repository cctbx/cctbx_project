from __future__ import division
from libtbx import adopt_init_args
from libtbx.str_utils import format_value
from cctbx import maptbx
from cctbx.array_family import flex
from mmtbx.rotamer.rotamer_eval import RotamerEval
from cctbx import miller
import mmtbx.monomer_library
from libtbx.utils import user_plus_sys_time
import mmtbx.refinement.real_space.individual_sites

time_target_simple = 0
time_compute_map   = 0
time_map_cc        = 0
time_update        = 0
time_states        = 0

def show_time(external):
  total = 0
  print "Detailed:"
  print "  time_target_simple: %6.3f" % time_target_simple
  print "  time_compute_map  : %6.3f" % time_compute_map
  print "  time_map_cc       : %6.3f" % time_map_cc
  print "  time_update       : %6.3f" % time_update
  print "  time_states       : %6.3f" % time_states
  total = time_target_simple+\
          time_compute_map+\
          time_map_cc+\
          time_update+\
          time_states
  for e in external:
    print e[0]%e[1]
    total += e[1]
  print "    sub-total       : %6.3f" % total

class states(object):
  def __init__(self, xray_structure, pdb_hierarchy):
    adopt_init_args(self, locals())
    self.counter = 0
    self.root = iotbx.pdb.hierarchy.root()
    self.sites_carts = []

  def add(self, sites_cart):
    global time_states
    timer = user_plus_sys_time()
    self.sites_carts.append(sites_cart)
    ph = self.pdb_hierarchy.deep_copy()
    xrs = self.xray_structure.replace_sites_cart(new_sites = sites_cart)
    ph.adopt_xray_structure(xrs)
    models = ph.models()
    md = models[0].detached_copy()
    md.id = str(self.counter)
    self.root.append_model(md)
    self.counter += 1
    time_states += timer.elapsed()

  def write(self, file_name):
    self.root.write_pdb_file(file_name = file_name)

def target_simple(target_map, sites_cart=None, sites_frac=None, unit_cell=None):
  global time_target_simple
  timer = user_plus_sys_time()
  assert [sites_cart, sites_frac].count(None) == 1
  if(sites_frac is None):
    assert unit_cell is not None
    sites_frac = unit_cell.fractionalize(sites_cart)
  result = 0
  for site_frac in sites_frac:
    result += target_map.eight_point_interpolation(site_frac)
  scale = 1
  if(sites_cart.size() > 0): scale = 1./sites_cart.size()
  time_target_simple += timer.elapsed()
  return result*scale

class evaluator(object):
  def __init__(self, sites_cart, target_map, unit_cell, target):
    adopt_init_args(self, locals())
    if(target is not None):
      self.target = -9999
    else:
      self.target = target_simple(
        target_map = self.target_map,
        sites_cart = self.sites_cart,
        unit_cell  = self.unit_cell)

  def evaluate(self, sites_cart):
    t = target_simple(
      target_map = self.target_map,
      sites_cart = sites_cart,
      unit_cell  = self.unit_cell)
    #print "  ", t
    if(t >= self.target):
      self.sites_cart = sites_cart
      self.target = t

class rsr_residue(object):
  def __init__(self,
               pdb_hierarchy_residue,
               selection_sidechain,
               selection_backbone,
               selection_all,
               map_cc_sidechain=None,
               map_cc_backbone=None,
               map_cc_all=None,
               map_value_sidechain=None,
               map_value_backbone=None,
               distance_to_closest_rotamer=None,
               rotamer_status=None):
    adopt_init_args(self, locals())

def rotamer_fit(residue, target_map, mon_lib_srv, unit_cell, rotamer_manager,
                xray_structure, residue_selection):
  import mmtbx.refinement.real_space.fit_residue
  mmtbx.refinement.real_space.fit_residue.manager(
    target_map  = target_map,
    mon_lib_srv = mon_lib_srv,
    unit_cell   = unit_cell,
    residue     = residue)
  return residue.atoms().extract_xyz()

def get_rotamer_iterator(mon_lib_srv, residue):
  get_class = iotbx.pdb.common_residue_names_get_class
  rotamer_iterator = None
  if(get_class(residue.resname) == "common_amino_acid"):
    rotamer_iterator = mon_lib_srv.rotamer_iterator(
      fine_sampling = True,
      comp_id       = residue.resname,
      atom_names    = residue.atoms().extract_name(),
      sites_cart    = residue.atoms().extract_xyz())
    if(rotamer_iterator is None or
       rotamer_iterator.problem_message is not None or
       rotamer_iterator.rotamer_info is None):
      rotamer_iterator = None
  return rotamer_iterator

class monitor(object):
  def __init__(self,
               pdb_hierarchy,
               xray_structure,
               target_map,
               geometry_restraints_manager,
               xray_structure_reference = None):
    adopt_init_args(self, locals())
    self.unit_cell = self.xray_structure.unit_cell()
    self.xray_structure_start = xray_structure.deep_copy_scatterers()
    self.states_collector = states(
      pdb_hierarchy  = self.pdb_hierarchy,
      xray_structure = self.xray_structure)
    self.cc = None
    self.rmsd_b = None
    self.rmsd_a = None
    self.dist_from_start = None
    self.dist_from_ref = None
    self.number_of_rotamer_outliers = 0
    self.rotamer_manager = RotamerEval()
    self.residues = self.get_residues()
    self.set_globals()
    self.mon_lib_srv = mmtbx.monomer_library.server.server()
    self.set_rsr_residue_attributes()

  def compute_map(self, xray_structure):
    global time_compute_map
    timer = user_plus_sys_time()
    mc = self.target_map.miller_array.structure_factors_from_scatterers(
      xray_structure = xray_structure).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.target_map.crystal_gridding,
      fourier_coefficients = mc)
    fft_map.apply_sigma_scaling()
    time_compute_map += timer.elapsed()
    return fft_map.real_map_unpadded()

  def map_cc(self, map, sites_cart = None):
    global time_map_cc
    timer = user_plus_sys_time()
    result = None
    if(sites_cart is not None):
      sel = maptbx.grid_indices_around_sites(
        unit_cell  = self.unit_cell,
        fft_n_real = map.focus(),
        fft_m_real = map.all(),
        sites_cart = sites_cart,
        site_radii = flex.double(sites_cart.size(), 2))
      result = flex.linear_correlation(
        x=map.select(sel).as_1d(),
        y=self.target_map.data.select(sel).as_1d()).coefficient()
    else:
      result = flex.linear_correlation(
        x=map.as_1d(),
        y=self.target_map.data.as_1d()).coefficient()
    time_map_cc += timer.elapsed()
    return result

  def show_residues(self):
    print
    print "resid    CC(sc) CC(bb) CC(sb)   Rotamer  Dist. to nearest rotamer"
    fmt="%s %s %6.3f %6.3f %6.3f %9s %s"
    for r in self.residues:
      ms=r.map_cc_sidechain
      mb=r.map_cc_backbone
      ma=r.map_cc_all
      if(ms<-1): ms = -1
      if(mb<-1): mb = -1
      if(ma<-1): ma = -1
      print fmt % (
        r.pdb_hierarchy_residue.resname,
        r.pdb_hierarchy_residue.resseq,
        ms,
        mb,
        ma,
        str(r.rotamer_status),
        format_value("%6.3f",r.distance_to_closest_rotamer))

  def set_rsr_residue_attributes(self):
    small = -1.e9
    get_class = iotbx.pdb.common_residue_names_get_class
    unit_cell = self.xray_structure.unit_cell()
    current_map = self.compute_map(xray_structure = self.xray_structure)
    sites_cart = self.xray_structure.sites_cart()
    for r in self.residues:
      sca = sites_cart.select(r.selection_all)
      scs = sites_cart.select(r.selection_sidechain)
      scb = sites_cart.select(r.selection_backbone)
      r.rotamer_status = self.rotamer_manager.evaluate_residue(
        residue=r.pdb_hierarchy_residue)
      if(r.rotamer_status != "OUTLIER"):
        r.map_cc_all       = self.map_cc(sites_cart = sca, map = current_map)
        r.map_cc_sidechain = self.map_cc(sites_cart = scs, map = current_map)
        r.map_cc_backbone  = self.map_cc(sites_cart = scb, map = current_map)
        r.map_value_sidechain = target_simple(target_map=current_map,
          sites_cart=scs, unit_cell=self.unit_cell)
        r.map_value_backbone = target_simple(target_map=current_map,
          sites_cart=scb, unit_cell=self.unit_cell)
      else:
        r.map_cc_all          = small
        r.map_cc_sidechain    = small
        r.map_cc_backbone     = small
        r.map_value_sidechain = small
        r.map_value_backbone  = small
      rotamer_iterator = get_rotamer_iterator(
        mon_lib_srv = self.mon_lib_srv,
        residue     = r.pdb_hierarchy_residue)
      if(rotamer_iterator is not None):
        dist_min = 1.e+9
        for ro, rotamer_sites_cart in rotamer_iterator:
          d=flex.mean(flex.sqrt((sites_cart.select(r.selection_all) - rotamer_sites_cart).dot()))
          if(d < dist_min):
            dist_min = d
        r.distance_to_closest_rotamer = dist_min

  def get_residues(self):
    residue_groups = self.pdb_hierarchy.residue_groups()
    backbone_atoms = ["N","CA","C","O","CB"]
    selections_all = []
    result = []
    for residue_group in residue_groups:
      for conformer in residue_group.conformers():
        for residue in conformer.residues():
          residue_i_seqs_backbone  = flex.size_t()
          residue_i_seqs_sidechain = flex.size_t()
          residue_i_seqs_all = flex.size_t()
          for atom in residue.atoms():
            an = atom.name.strip()
            bb = an in backbone_atoms
            residue_i_seqs_all.append(atom.i_seq)
            if(bb): residue_i_seqs_backbone.append(atom.i_seq)
            else:   residue_i_seqs_sidechain.append(atom.i_seq)
          selections_all.append(residue_i_seqs_all)
          result.append(rsr_residue(
            pdb_hierarchy_residue = residue,
            selection_sidechain   = residue_i_seqs_sidechain,
            selection_backbone    = residue_i_seqs_backbone,
            selection_all         = residue_i_seqs_all))
    sa = flex.size_t()
    for s in selections_all: sa.extend(s)
    assert self.xray_structure.select(sa).scatterers().size() == \
           self.xray_structure.scatterers().size()
    return result

  def set_globals(self):
    self.cc = self.map_cc(map=self.compute_map(
      xray_structure = self.xray_structure))
    es = self.geometry_restraints_manager.energies_sites(
      sites_cart = self.xray_structure.sites_cart())
    self.rmsd_a = es.angle_deviations()[2]
    self.rmsd_b = es.bond_deviations()[2]
    self.dist_from_start = flex.mean(self.xray_structure_start.distances(
      other = self.xray_structure))
    if(self.xray_structure_reference is not None):
      self.dist_from_ref = flex.mean(self.xray_structure_reference.distances(
        other = self.xray_structure))
    self.number_of_rotamer_outliers = 0
    for r in self.residues:
      rotamer_status = self.rotamer_manager.evaluate_residue(
        residue=r.pdb_hierarchy_residue)
      if(rotamer_status == "OUTLIER"):
        self.number_of_rotamer_outliers += 1

  def update(self, xray_structure, accept_any=False):
    global time_update
    timer = user_plus_sys_time()
    unit_cell = xray_structure.unit_cell()
    current_map = self.compute_map(xray_structure = xray_structure)
    sites_cart  = xray_structure.sites_cart()
    sites_cart_ = self.xray_structure.sites_cart()
    for r in self.residues:
      sca = sites_cart.select(r.selection_all)
      scs = sites_cart.select(r.selection_sidechain)
      scb = sites_cart.select(r.selection_backbone)
      map_cc_all       = self.map_cc(sites_cart = sca, map = current_map)
      map_cc_sidechain = self.map_cc(sites_cart = scs, map = current_map)
      map_cc_backbone  = self.map_cc(sites_cart = scb, map = current_map)
      map_value_sidechain = target_simple(target_map=current_map,
        sites_cart=scs, unit_cell=self.unit_cell)
      map_value_backbone = target_simple(target_map=current_map,
        sites_cart=scb, unit_cell=self.unit_cell)
      flag = map_cc_all      > r.map_cc_all and \
             map_cc_backbone > r.map_cc_backbone and \
             map_cc_backbone > map_cc_sidechain
      if(r.map_value_backbone > r.map_value_sidechain):
        if(map_value_backbone < map_value_sidechain):
          flag = False
      if(accept_any): flag=True
      if(flag):
        residue_sites_cart_new = sites_cart.select(r.selection_all)
        sites_cart_.set_selected(r.selection_all, residue_sites_cart_new)
        r.pdb_hierarchy_residue.atoms().set_xyz(residue_sites_cart_new)
        rotamer_status = self.rotamer_manager.evaluate_residue(
          residue=r.pdb_hierarchy_residue)
        r.rotamer_status = rotamer_status
    self.xray_structure= self.xray_structure.replace_sites_cart(sites_cart_)
    self.set_globals()
    self.set_rsr_residue_attributes()
    self.states_collector.add(sites_cart = sites_cart_)
    time_update += timer.elapsed()

  def show(self, suffix=""):
    if(self.dist_from_ref is None):
      f="cc: %6.4f rmsd_b: %6.4f rmsd_a: %5.2f d(start): %6.3f rota: %3d"
      print f%(self.cc,self.rmsd_b,self.rmsd_a,self.dist_from_start,
               self.number_of_rotamer_outliers),suffix
    else:
      f="cc: %6.4f rmsd_b: %6.4f rmsd_a: %5.2f d(start): %6.3f d(ref): %6.3f"
      print f%(self.cc,self.rmsd_b,self.rmsd_a,self.dist_from_start,
        self.dist_from_ref),suffix

def run(target_map,
        pdb_hierarchy,
        xray_structure,
        geometry_restraints_manager,
        xray_structure_reference = None,
        rms_bonds_limit  = 0.03,
        rms_angles_limit = 3.0,
        max_iterations   = 500,
        macro_cycles     = 20,
        minimization     = True,
        expload          = True,
        rotamer_search   = True,
        verbose          = True):
  sel = flex.bool(xray_structure.scatterers().size(), True)
  d_min = target_map.miller_array.d_min()
  #
  #geometry_restraints_manager.geometry.remove_dihedrals_in_place(sel)
  #
  rsr_simple_refiner = mmtbx.refinement.real_space.individual_sites.simple(
    target_map                  = target_map.data,
    selection                   = sel,
    real_space_gradients_delta  = d_min/4,
    max_iterations              = max_iterations,
    geometry_restraints_manager = geometry_restraints_manager.geometry)

  #
  if(xray_structure_reference is not None):
    xray_structure_reference = xray_structure_reference.deep_copy_scatterers()
  monitor_object = monitor(
    pdb_hierarchy = pdb_hierarchy,
    xray_structure = xray_structure.deep_copy_scatterers(),
    target_map = target_map,
    geometry_restraints_manager = geometry_restraints_manager.geometry,
    xray_structure_reference = xray_structure_reference)
  if(verbose):
    monitor_object.show(suffix="start")
    monitor_object.show_residues()
  #
  tmp = monitor_object.xray_structure.deep_copy_scatterers()
  weight_d = 50
  weight_s = 50
  #
  weights = flex.double()
  optimize_weight = True
  ####
  #sites_cart = tmp.sites_cart()
  #for i_r, r in enumerate(monitor_object.residues):
  #  print i_r, len(monitor_object.residues)
  #  sites_cart_ = rotamer_fit(
  #    residue     = r.pdb_hierarchy_residue,
  #    target_map  = target_map.data,
  #    mon_lib_srv = monitor_object.mon_lib_srv,
  #    xray_structure = tmp,
  #    residue_selection = r.selection_all,
  #    unit_cell   = tmp.unit_cell(),
  #    rotamer_manager = monitor_object.rotamer_manager)
  #  sites_cart.set_selected(r.selection_all, sites_cart_)
  #tmp = tmp.replace_sites_cart(sites_cart)
  #monitor_object.update(xray_structure=tmp, accept_any=True)
  #monitor_object.show(suffix="start: initial rotamer fixing")
  #
  #pdb_hierarchy.adopt_xray_structure(tmp)
  #pdb_hierarchy.write_pdb_file(file_name = "ZZZ.pdb")
  #
  #STOP()
  ####
  geometry_restraints_manager.geometry.generic_restraints_manager.reference_manager.\
    remove_chi_angle_restraints(pdb_hierarchy=pdb_hierarchy)
  geometry_restraints_manager.geometry.generic_restraints_manager.reference_manager.\
    add_torsion_restraints(
      pdb_hierarchy   = pdb_hierarchy,
      sites_cart      = tmp.sites_cart(),
      chi_angles_only = True,
      sigma           = 10.0)
  ####
  for i in range(macro_cycles):
    tmp_dc = tmp.deep_copy_scatterers()
    if(expload and i>1 and i%2==0):
      tmp.shake_sites_in_place(mean_distance=3) # reverse back if refinement failed
      #
      #from mmtbx.dynamics import cartesian_dynamics
      #grad_calc = cartesian_dynamics.gradients_calculator_geometry_restraints(
      #  restraints_manager = geometry_restraints_manager.geometry)
      #cartesian_dynamics.run(
      #  xray_structure       = tmp,
      #  gradients_calculator = grad_calc,
      #  temperature          = 3000,
      #  n_steps              = 100000,
      #  time_step            = 0.0005,
      #  stop_cm_motion       = True,
      #  stop_at_diff         = 5.0)
      #
      #from mmtbx.dynamics import simulated_annealing as sa
      #params = sa.master_params().extract()
      #params.start_temperature=10000
      #params.final_temperature=0
      #params.cool_rate = 1000
      #params.number_of_steps = 500
      #params.update_grads_shift = 0.
      #sa.run(
      #  params             = params,
      #  #fmodel             = fmodel,
      #  xray_structure     = tmp,
      #  real_space         = True,
      #  target_map         = target_map.data,
      #  restraints_manager = geometry_restraints_manager,
      #  wx                 = weight_s,
      #  wc                 = 1.,
      #  verbose            = True)
      #
    if(minimization):
      target_type = "simple"
      refined = mmtbx.refinement.real_space.individual_sites.refinery(
        refiner          = rsr_simple_refiner,
        optimize_weight  = optimize_weight,
        xray_structure   = tmp,
        start_trial_weight_value = weight_s,
        rms_bonds_limit  = rms_bonds_limit,
        rms_angles_limit = rms_angles_limit)
      if(refined.sites_cart_result is not None):
        tmp = tmp.replace_sites_cart(refined.sites_cart_result)
        weight_s = refined.weight_final
        monitor_object.update(xray_structure=tmp)#, accept_any=True)
        if(verbose):
          monitor_object.show(suffix=" weight: %s"%str(weight_s))
        tmp = monitor_object.xray_structure.deep_copy_scatterers()
        weights.append(refined.weight_final)
        if(weights.size() == 2):
          weight_s = flex.mean(weights)
          optimize_weight = False
      else:
        tmp = tmp_dc
        print "Refinement failed."
        #
    #if(rotamer_search or ((not expload or (expload and minimization)) and i>macro_cycles/2)):
    if(not (i>1 and i%2==0) and (rotamer_search or (not expload or (expload and minimization)))):
      sites_cart = tmp.sites_cart()
      for r in monitor_object.residues:
        sites_cart_ = rotamer_fit(
          residue     = r.pdb_hierarchy_residue,
          target_map  = target_map.data,
          mon_lib_srv = monitor_object.mon_lib_srv,
          unit_cell   = tmp.unit_cell(),
          xray_structure = tmp,
          residue_selection = r.selection_all,
          rotamer_manager = monitor_object.rotamer_manager)
        sites_cart.set_selected(r.selection_all, sites_cart_)
      tmp = tmp.replace_sites_cart(sites_cart)
      monitor_object.update(xray_structure=tmp, accept_any=True)
      # add reference restraints
      geometry_restraints_manager.geometry.generic_restraints_manager.reference_manager.\
        remove_chi_angle_restraints(pdb_hierarchy=pdb_hierarchy)
      geometry_restraints_manager.geometry.generic_restraints_manager.reference_manager.\
        add_torsion_restraints(
          pdb_hierarchy   = pdb_hierarchy,
          sites_cart      = tmp.sites_cart(),
          chi_angles_only = True,
          sigma           = 2.0)
      #
      if(verbose):
        monitor_object.show(suffix=" weight: %s"%str(None))
      tmp = monitor_object.xray_structure.deep_copy_scatterers()
        #

  #
  if(minimization):
    refined = refinery(
      optimize_weight  = True,
      refiner          = rsr_simple_refiner,
      xray_structure   = tmp,
      start_trial_weight_value = 50,
      rms_bonds_limit  = 0.02,
      rms_angles_limit = 3.0)
    if(verbose):
      print "FINAL:", refined.rms_bonds_final,refined.rms_angles_final
    if(refined.sites_cart_result is not None):
      tmp = tmp.replace_sites_cart(refined.sites_cart_result)
      weight_s = refined.weight_final
      monitor_object.update(xray_structure=tmp, accept_any=True) # XXX ???
      if(verbose): monitor_object.show(suffix=" weight: %s"%str(weight_s))
  #
  if(verbose): monitor_object.show_residues()
  monitor_object.states_collector.write(file_name = "all.pdb")
  return monitor_object.xray_structure
