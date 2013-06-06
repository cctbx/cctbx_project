from __future__ import division
from libtbx import adopt_init_args
from cctbx.array_family import flex
import mmtbx.monomer_library
from libtbx.utils import user_plus_sys_time
import mmtbx.refinement.real_space.individual_sites
from mmtbx import monomer_library
import mmtbx.refinement.real_space.fit_residues

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

def morph(structure_monitor):
  sm = structure_monitor
  from solve_resolve.resolve_python import resolve_in_memory
  import iotbx.pdb
  input_text="""morph  MEMORY
morph_main
rad_morph 6
      """
  sm.pdb_hierarchy.write_pdb_file(file_name="tmp_morphing_in.pdb")
  pdb_inp = iotbx.pdb.input(file_name = "tmp_morphing_in.pdb")
  cmn=resolve_in_memory.run(
    input_text=input_text,
    map_coeffs=sm.target_map_object.miller_array,
    model_pdb_inp=pdb_inp)
  morphed_model_pdb=cmn.atom_db.pdb_out_as_string
  of = open("tmp_morphed.pdb","w")
  print >> of, morphed_model_pdb
  of.close()
  pdb_inp = iotbx.pdb.input(file_name = "tmp_morphed.pdb")
  pdb_hierarchy = pdb_inp.construct_hierarchy()
  pdb_hierarchy.atoms().reset_i_seq()
  xrs = pdb_inp.xray_structure_simple()
  sm.update(xray_structure=xrs, accept_as_is=True)
  sm.show(prefix="morph:")
  return sm.xray_structure.deep_copy_scatterers()

def simulated_annealing(structure_monitor, tmp):
  sm = structure_monitor
  from mmtbx.dynamics import simulated_annealing as sa
  params = sa.master_params().extract()
  params.start_temperature=5000
  params.final_temperature=0
  params.cool_rate = 300
  params.number_of_steps = 100
  params.update_grads_shift = 0.
  sa.run(
    params             = params,
    xray_structure     = tmp,
    real_space         = True,
    target_map         = target_map.data,
    restraints_manager = geometry_restraints_manager,
    wx                 = weight_s,
    wc                 = 1.,
    verbose            = True)
  sm.update(xray_structure=tmp, accept_as_is=True)
  sm.show(prefix="SA:")

def refine_residue(structure_monitor, mon_lib_srv, geometry_restraints_manager):
  sm = structure_monitor
  result = mmtbx.refinement.real_space.fit_residues.manager(
    structure_monitor = sm,
    mon_lib_srv       = mon_lib_srv)
  sm.show(prefix="rota_s:")
  tmp = sm.xray_structure.deep_copy_scatterers()
  def switch_to_exact_rotamers(tmp, sm):
    xrs = tmp.deep_copy_scatterers()
    ph = sm.pdb_hierarchy.deep_copy()
    ph.atoms().reset_i_seq()
    xrs = mmtbx.utils.max_distant_rotomer(
      xray_structure = xrs,
      pdb_hierarchy  = ph,
      exact_match    = True)
    return ph, xrs
  ph, xrs = switch_to_exact_rotamers(tmp=tmp, sm=sm)
  geometry_restraints_manager.geometry.generic_restraints_manager.reference_manager.\
    remove_chi_angle_restraints(pdb_hierarchy=sm.pdb_hierarchy)
  geometry_restraints_manager.geometry.generic_restraints_manager.reference_manager.\
    add_torsion_restraints(
      pdb_hierarchy   = ph, #sm.pdb_hierarchy,
      sites_cart      = xrs.sites_cart(),#tmp.sites_cart(),
      chi_angles_only = True,
      sigma           = 5.0)
  return tmp

def run(target_map,
        pdb_hierarchy,
        xray_structure,
        geometry_restraints_manager,
        xray_structure_reference = None,
        rms_bonds_limit         = 0.02,
        rms_angles_limit        = 2.0,
        max_iterations          = 500,
        macro_cycles            = 20,
        minimization            = True,
        rotamer_search          = True,
        verbose                 = True,
        run_morphing            = True,
        run_simulated_annealing = False,
        expload_size            = 1.0):
  sel = flex.bool(xray_structure.scatterers().size(), True)
  d_min = target_map.miller_array.d_min()
  mon_lib_srv = monomer_library.server.server()
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
  sm = mmtbx.refinement.real_space.structure_monitor(
    pdb_hierarchy               = pdb_hierarchy,
    xray_structure              = xray_structure.deep_copy_scatterers(),
    target_map_object           = target_map,
    geometry_restraints_manager = geometry_restraints_manager.geometry)
  if(verbose):
    sm.show_residues()
    sm.show(prefix="start:")
  #
  tmp = sm.xray_structure.deep_copy_scatterers()
  weight_d = 50
  weight_s = 50
  #
  weights = flex.double()
  optimize_weight = True
  ####
  expload_step = expload_size/macro_cycles

  for i in range(macro_cycles):
    tmp_dc = tmp.deep_copy_scatterers()
    if(i==0):
      tmp = refine_residue(structure_monitor=sm, mon_lib_srv=mon_lib_srv,
        geometry_restraints_manager=geometry_restraints_manager)
    if(i==0 and run_morphing):
      tmp = morph(structure_monitor = sm)
    if(run_simulated_annealing):
      simulated_annealing(structure_monitor = sm, tmp=tmp)
    if(expload_size is not None):
      tmp.shake_sites_in_place(mean_distance=expload_size)
      expload_size = max(expload_size-expload_step, 0)

    if(minimization):
      refined = mmtbx.refinement.real_space.individual_sites.refinery(
        refiner                  = rsr_simple_refiner,
        optimize_weight          = True,#optimize_weight,
        xray_structure           = tmp,
        start_trial_weight_value = weight_s,
        rms_bonds_limit          = rms_bonds_limit,
        rms_angles_limit         = rms_angles_limit)
      if(refined.sites_cart_result is not None):
        tmp = tmp.replace_sites_cart(refined.sites_cart_result)
        weight_s = refined.weight_final
        sm.update(xray_structure=tmp, accept_as_is=True) #XXX
        if(verbose):
          sm.show(prefix="mc %s:"%str(weight_s))
        tmp = sm.xray_structure.deep_copy_scatterers()
        weights.append(refined.weight_final)
        if(weights.size() == 2):
          weight_s = flex.mean(weights)
          optimize_weight = False
      else:
        tmp = tmp_dc
        print "Refinement failed."
      #
      tmp = refine_residue(structure_monitor=sm, mon_lib_srv=mon_lib_srv,
        geometry_restraints_manager=geometry_restraints_manager)
  #
  if(minimization):
    refined = mmtbx.refinement.real_space.individual_sites.refinery(
      optimize_weight          = True,
      refiner                  = rsr_simple_refiner,
      xray_structure           = tmp,
      start_trial_weight_value = 50,
      rms_bonds_limit          = 0.015,
      rms_angles_limit         = 1.5)
    if(verbose):
      print "FINAL:", refined.rms_bonds_final,refined.rms_angles_final
    if(refined.sites_cart_result is not None):
      tmp = tmp.replace_sites_cart(refined.sites_cart_result)
      weight_s = refined.weight_final
      sm.update(xray_structure=tmp, accept_as_is=True) # XXX ???
      if(verbose): sm.show(prefix="%s"%str(weight_s))
  #
  if(verbose): sm.show_residues()
  sm.states_collector.write(file_name = "all.pdb")
  return sm.xray_structure
