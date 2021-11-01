from __future__ import absolute_import, division, print_function
from libtbx import adopt_init_args
import scitbx.lbfgs
from cctbx import maptbx
from cctbx.array_family import flex
from mmtbx import utils
from libtbx.test_utils import approx_equal
from cctbx import crystal
import mmtbx.refinement.minimization_ncs_constraints
import mmtbx.model
from cctbx.maptbx import minimization
from six.moves import range

class easy(object):
  """
  Simplest interface to most automated and fast real-space refinement.
  To keep it simple not all parameters are exposed.
  """
  def __init__(
        self,
        map_data,
        xray_structure,
        pdb_hierarchy,
        geometry_restraints_manager,
        gradients_method="fd",
        selection=None,
        selection_real_space=None,
        rms_bonds_limit=0.015,
        rms_angles_limit=2.0,
        max_iterations=150,
        w = None,
        states_accumulator=None,
        log=None):
    assert gradients_method in ["fd", "linear", "quadratic", "tricubic"]
    adopt_init_args(self, locals())
    es = geometry_restraints_manager.geometry.energies_sites(
      sites_cart = xray_structure.sites_cart())
    self.rmsd_angles_start = es.angle_deviations()[2]
    self.rmsd_bonds_start  = es.bond_deviations()[2]
    self.w = w
    if(self.max_iterations > 0):
      if(self.w is None):
        import mmtbx.refinement.real_space.weight
        self.weight = mmtbx.refinement.real_space.weight.run(
          map_data                    = map_data,
          xray_structure              = self.xray_structure,
          pdb_hierarchy               = self.pdb_hierarchy,
          geometry_restraints_manager = geometry_restraints_manager,
          rms_bonds_limit             = rms_bonds_limit,
          rms_angles_limit            = rms_angles_limit,
          gradients_method            = gradients_method)
        self.w = self.weight.weight
      if(selection is None):
        selection = flex.bool(self.xray_structure.scatterers().size(), True)
      refine_object = simple(
        target_map                  = map_data,
        selection                   = selection,
        selection_real_space        = selection_real_space,
        max_iterations              = max_iterations,
        geometry_restraints_manager = geometry_restraints_manager.geometry,
        states_accumulator          = states_accumulator,
        gradients_method            = gradients_method)
      refine_object.refine(
        weight = self.w, xray_structure = self.xray_structure)
      self.rmsd_bonds_final, self.rmsd_angles_final = refine_object.rmsds()
      self.xray_structure=self.xray_structure.replace_sites_cart(
        new_sites=refine_object.sites_cart(), selection=None)
      self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

class simple(object):
  def __init__(
        self,
        target_map,
        selection,
        geometry_restraints_manager,
        gradients_method="fd",
        real_space_gradients_delta=1./4,
        selection_real_space=None,
        max_iterations=150,
        states_accumulator=None,
        ncs_groups=None):
    adopt_init_args(self, locals())
    self.lbfgs_core_params = scitbx.lbfgs.core_parameters(
      stpmax = 25.0)
    self.lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    self.lbfgs_exception_handling_params = scitbx.lbfgs.\
      exception_handling_parameters(
        ignore_line_search_failed_step_at_lower_bound = True,
        ignore_line_search_failed_step_at_upper_bound = True,
        ignore_line_search_failed_maxfev              = True)
    self.refined = None
    self.crystal_symmetry = None
    self.site_symmetry_table = None

  def refine(self, weight, xray_structure):
    if(self.ncs_groups is None):
      self.crystal_symmetry = xray_structure.crystal_symmetry()
      self.site_symmetry_table = xray_structure.site_symmetry_table()
      self.refined = maptbx.real_space_refinement_simple.lbfgs(
        selection_variable              = self.selection,
        gradients_method                = self.gradients_method,
        selection_variable_real_space   = self.selection_real_space,
        sites_cart                      = xray_structure.sites_cart(),
        density_map                     = self.target_map,
        geometry_restraints_manager     = self.geometry_restraints_manager,
        real_space_target_weight        = weight,
        real_space_gradients_delta      = self.real_space_gradients_delta,
        lbfgs_core_params               = self.lbfgs_core_params,
        lbfgs_termination_params        = self.lbfgs_termination_params,
        lbfgs_exception_handling_params = self.lbfgs_exception_handling_params,
        states_collector                = self.states_accumulator)
    else:
      refine_selection = flex.bool(xray_structure.scatterers().size(),
        True).iselection()
      self.refined = mmtbx.refinement.minimization_ncs_constraints.\
        target_function_and_grads_real_space(
          map_data                   = self.target_map,
          xray_structure             = xray_structure.deep_copy_scatterers(),
          ncs_restraints_group_list  = self.ncs_groups,
          refine_selection           = refine_selection,
          real_space_gradients_delta = self.real_space_gradients_delta,
          restraints_manager         = self.geometry_restraints_manager,
          data_weight                = weight,
          refine_sites               = True)
      minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
        target_and_grads_object      = self.refined,
        xray_structure               = xray_structure.deep_copy_scatterers(),
        ncs_restraints_group_list    = self.ncs_groups,
        refine_selection             = refine_selection,
        finite_grad_differences_test = False,
        max_iterations               = self.max_iterations,
        refine_sites                 = True)

  def sites_cart(self):
    assert self.refined is not None
    if(self.ncs_groups is None):
      sites_cart = self.refined.sites_cart
    else:
      sites_cart = self.refined.xray_structure.sites_cart()
      return sites_cart #XXX TRAP
    if(self.selection):
      sites_cart.set_selected(self.selection, self.refined.sites_cart_variable)
    special_position_indices = self.site_symmetry_table.special_position_indices()
    if(special_position_indices.size()>0):
      for i_seq in special_position_indices:
        sites_cart[i_seq] = crystal.correct_special_position(
          crystal_symmetry = self.crystal_symmetry,
          special_op       = self.site_symmetry_table.get(i_seq).special_op(),
          site_cart        = sites_cart[i_seq],
          site_label       = None,
          tolerance        = 1)
    return sites_cart

  def rmsds(self):
    b,a = None,None
    es = self.geometry_restraints_manager.energies_sites(
      sites_cart = self.sites_cart())
    a = es.angle_deviations()[2]
    b = es.bond_deviations()[2]
    return b,a

class diff_map(object):
  def __init__(self,
               miller_array,
               crystal_gridding,
               map_target,
               geometry_restraints_manager,
               restraints_target_weight = 1,
               max_iterations = 500,
               min_iterations = 500):
    adopt_init_args(self, locals())
    self.step = miller_array.d_min()/4.
    self.refined = None

  def refine(self, weight, sites_cart=None, xray_structure=None):
    assert xray_structure is not None and [sites_cart,xray_structure].count(None)==1
    self.refined = minimization.run(
      xray_structure              = xray_structure,
      miller_array                = self.miller_array,
      crystal_gridding            = self.crystal_gridding,
      map_target                  = self.map_target,
      max_iterations              = self.max_iterations,
      min_iterations              = self.min_iterations,
      step                        = self.step,
      real_space_target_weight    = weight,
      restraints_target_weight    = self.restraints_target_weight,
      geometry_restraints_manager = self.geometry_restraints_manager,
      target_type                 = "diffmap")

  def sites_cart(self):
    assert self.refined is not None
    return self.refined.xray_structure.sites_cart()

class refinery(object):
  def __init__(self,
               refiner,
               xray_structure,
               start_trial_weight_value = 50.,
               weight_sample_rate = 10,
               rms_bonds_limit = 0.03,
               rms_angles_limit = 3.0,
               optimize_weight = True):
    self.rms_angles_start = None
    self.rms_bonds_start = None
    self.refiner = refiner
    self.weight_start=start_trial_weight_value
    sites_cart_start = xray_structure.sites_cart()
    self.rms_bonds_start, self.rms_angles_start  = \
      self.rmsds(sites_cart=xray_structure.sites_cart())
    self.weight_sample_rate = weight_sample_rate
    # results
    self.weight_final = None
    self.sites_cart_result = None
    self.rms_bonds_final,self.rms_angles_final = None,None
    #
    pool = {}
    bonds = flex.double()
    angles = flex.double()
    weights = flex.double()
    #
    weight = start_trial_weight_value
    weight_last = weight
    self.adjust_weight_sample_rate(weight=weight)
    if(optimize_weight):
      while True:
        self.rmsds(sites_cart=sites_cart_start) # DUMMY
        self.adjust_weight_sample_rate(weight=weight_last)
        tmp = xray_structure.deep_copy_scatterers()
        #tmp.shake_sites_in_place(
        #  rms_difference = None,
        #  mean_distance  = 0.5)
        refiner.refine(
          xray_structure = tmp,#xray_structure.deep_copy_scatterers(), # XXX
          weight     = weight)
        sites_cart_result = refiner.sites_cart()
        bd, ad = self.rmsds(sites_cart=sites_cart_result)
        bonds.append(bd)
        angles.append(ad)
        weights.append(weight)
        pool.setdefault(weight,[]).append([sites_cart_result.deep_copy(),bd,ad])
        if(refiner.geometry_restraints_manager is None): break
        weight_last = weight
        if(ad>rms_angles_limit or bd > rms_bonds_limit):
          weight -= self.weight_sample_rate
        else:
          weight += self.weight_sample_rate
        if(weight<0 or abs(weight)<1.e-6):
          self.adjust_weight_sample_rate(weight=weight)
          weight = weight_last
          weight -= self.weight_sample_rate
        #print ">>> ", "%8.4f %8.4f"%(weight, weight_last), "%6.4f %5.2f"%(bd, ad),\
        #  self.weight_sample_rate, "  f (start/final):", refiner.refined.f_start, refiner.refined.f_final
        if((weight<0 or weight>1000) or weight in weights): break
        l = bonds.size()-1
        if(bonds.size()>5 and
           abs(bonds[l]-bonds[l-1])<0.0005 and
           abs(bonds[l]-bonds[l-2])<0.0005 and
           abs(bonds[l]-bonds[l-3])<0.0005 and
           abs(bonds[l]-bonds[l-4])<0.0005 and
           abs(bonds[l]-bonds[l-5])<0.0005): break
    else:
      refiner.refine(
        xray_structure = xray_structure.deep_copy_scatterers(), # XXX
        weight     = weight)
      sites_cart_result = refiner.sites_cart()
    # select results
    if(optimize_weight):
      delta = bonds-rms_bonds_limit
      ind = (delta == flex.max_default(delta.select(delta<=0),
        flex.min(delta))).iselection()[0]
      self.weight_final = weights[ind]
      self.sites_cart_result = pool[self.weight_final][0][0]
      self.rms_bonds_final,self.rms_angles_final = \
        self.rmsds(sites_cart=self.sites_cart_result)
      assert approx_equal(pool[self.weight_final][0][2], angles[ind])
      assert approx_equal(pool[self.weight_final][0][1], bonds[ind])
      assert approx_equal(self.rms_angles_final, angles[ind])
      assert approx_equal(self.rms_bonds_final, bonds[ind])
    else:
      self.weight_final = self.weight_start
      self.sites_cart_result = sites_cart_result

  def rmsds(self, sites_cart):
    b,a = None,None
    if(self.refiner.geometry_restraints_manager is not None):
      es = self.refiner.geometry_restraints_manager.energies_sites(
        sites_cart = sites_cart)
      a = es.angle_deviations()[2]
      b = es.bond_deviations()[2]
    return b,a

  def adjust_weight_sample_rate(self, weight):
    if(  weight <= 0.01 ): self.weight_sample_rate=0.001
    elif(weight <= 0.1  ): self.weight_sample_rate=0.01
    elif(weight <= 1.0  ): self.weight_sample_rate=0.1
    elif(weight <= 10.  ): self.weight_sample_rate=1.
    elif(weight <= 100. ): self.weight_sample_rate=10.
    elif(weight <= 1000.): self.weight_sample_rate=100.

class box_refinement_manager(object):
  def __init__(self,
               xray_structure,
               target_map,
               geometry_restraints_manager,
               gradients_method="fd",
               real_space_gradients_delta=1./4,
               max_iterations = 50,
               ncs_groups = None):
    self.xray_structure = xray_structure
    self.sites_cart = xray_structure.sites_cart()
    self.target_map = target_map
    self.geometry_restraints_manager = geometry_restraints_manager
    self.max_iterations=max_iterations
    self.real_space_gradients_delta = real_space_gradients_delta
    self.weight_optimal = None
    self.ncs_groups = ncs_groups
    self.gradients_method = gradients_method

  def update_xray_structure(self, new_xray_structure):
    self.xray_structure = new_xray_structure

  def update_target_map(self, new_target_map):
    self.target_map = new_target_map

  def refine(self,
             selection,
             optimize_weight = True,
             start_trial_weight_value = 50,
             selection_buffer_radius=5,
             box_cushion=2,
             rms_bonds_limit = 0.03,
             rms_angles_limit = 3.0):
    if(self.ncs_groups is None or len(self.ncs_groups)==0): # no NCS constraints
      sites_cart_moving = self.sites_cart
      selection_within = self.xray_structure.selection_within(
        radius    = selection_buffer_radius,
        selection = selection)
      sel = selection.select(selection_within)
      iselection = flex.size_t()
      for i, state in enumerate(selection):
        if state:
          iselection.append(i)
      box = utils.extract_box_around_model_and_map(
        xray_structure = self.xray_structure.select(selection=selection_within),
        map_data       = self.target_map,
        box_cushion    = box_cushion)
      new_unit_cell = box.xray_structure_box.unit_cell()
      geo_box = self.geometry_restraints_manager.select(selection_within)
      geo_box = geo_box.discard_symmetry(new_unit_cell=new_unit_cell)
      geo_box.shift_sites_cart(box.shift_cart) # disaster happens otherwise
      map_box = box.map_box
      sites_cart_box = box.xray_structure_box.sites_cart()

      rsr_simple_refiner = simple(
        target_map                  = map_box,
        selection                   = sel,
        real_space_gradients_delta  = self.real_space_gradients_delta,
        max_iterations              = self.max_iterations,
        geometry_restraints_manager = geo_box,
        gradients_method            = self.gradients_method)
      real_space_result = refinery(
        refiner                  = rsr_simple_refiner,
        xray_structure           = box.xray_structure_box,
        optimize_weight          = optimize_weight,
        start_trial_weight_value = start_trial_weight_value,
        rms_bonds_limit = rms_bonds_limit,
        rms_angles_limit = rms_angles_limit)
      self.weight_optimal = real_space_result.weight_final
      sites_cart_box_refined = real_space_result.sites_cart_result
      shift_back = [-box.shift_cart[i] for i in range(3)]
      sites_cart_box_refined_shifted_back = \
        sites_cart_box_refined + shift_back # Sure + not - ?
      sites_cart_refined = sites_cart_box_refined_shifted_back.select(sel)
      sites_cart_moving = sites_cart_moving.set_selected(
        iselection, sites_cart_refined)
      self.xray_structure.set_sites_cart(sites_cart_moving)
      self.sites_cart = self.xray_structure.sites_cart()
    else: # NCS constraints are present
      # select on xrs, grm, ncs_groups
      grm = self.geometry_restraints_manager.select(selection)
      xrs = self.xray_structure.select(selection)
      sel = flex.bool(xrs.scatterers().size(), True)
      size = self.xray_structure.scatterers().size()
      ncs_groups_ = self.ncs_groups.select(selection=selection)
      #
      rsr_simple_refiner = simple(
        target_map                  = self.target_map,
        selection                   = sel,
        real_space_gradients_delta  = self.real_space_gradients_delta,
        max_iterations              = self.max_iterations,
        geometry_restraints_manager = grm,
        ncs_groups                  = ncs_groups_)
      real_space_result = refinery(
        refiner                  = rsr_simple_refiner,
        xray_structure           = xrs,
        optimize_weight          = optimize_weight,
        start_trial_weight_value = start_trial_weight_value,
        rms_bonds_limit          = rms_bonds_limit,
        rms_angles_limit         = rms_angles_limit)
      self.weight_optimal = real_space_result.weight_final
#      # XXX undefined
#      self.xray_structure = None #XXX undefined
#      self.sites_cart = None     #XXX undefined


class minimize_wrapper_with_map():
  def __init__(self,
      model,
      target_map,
      refine_ncs_operators=False,
      number_of_cycles=1,
      cycles_to_converge=2,
      min_mode='simple_cycles',
      resolution=3.,
      log=None):

    # completely new way of doing this. using RSR macro-cycle
    # for test compatibility:
    print("Minimizing using reference map...", file=log)
    if model.ncs_constraints_present():
      print("  Minimizing... (NCS)", file=log)
    else:
      print("  Minimizing...", file=log)
    from phenix.refinement.macro_cycle_real_space import run as rsr_mc_run
    import scitbx.math
    from phenix.refinement import rsr
    rsr_master_params = rsr.master_params_str
    import iotbx.phil
    rsr_master_params = iotbx.phil.parse(rsr_master_params, process_includes=True)
    import mmtbx.idealized_aa_residues.rotamer_manager
    sin_cos_table = scitbx.math.sin_cos_table(n=10000)
    params = rsr_master_params.extract()
    params.pdb_interpretation = model._pdb_interpretation_params.pdb_interpretation
    params.refinement.run = "minimization_global+local_grid_search"
    params.refine_ncs_operators=False
    params.refinement.macro_cycles = number_of_cycles
    params.resolution = resolution
    rotamer_manager = mmtbx.idealized_aa_residues.rotamer_manager.load(
        rotamers = "favored")
    rigid_body_selections = [] # no RBR here
    from iotbx import map_model_manager
    from iotbx import map_manager
    mm = map_manager.map_manager(
      map_data                   = target_map,
      unit_cell_grid             = target_map.all(),
      wrapping                   = False,
      unit_cell_crystal_symmetry = model.crystal_symmetry())
    mmm = map_model_manager.map_model_manager(
      model       = model,
      map_manager = mm)
    res = rsr_mc_run(
      params                = params,
      map_model_manager     = mmm,
      log                   = log,
      rotamer_manager       = rotamer_manager,
      sin_cos_table         = sin_cos_table,
      rigid_body_selections = rigid_body_selections)

    model.set_sites_cart_from_hierarchy(res.model.get_hierarchy())
    res.structure_monitor.states_collector.write(file_name="rsr_all_states.pdb")
    return





    # end ===================================================
    # Very sophisticated implementation. Need to investigate before removal.

    from mmtbx.refinement.geometry_minimization import add_rotamer_restraints
    from mmtbx.refinement.minimization_monitor import minimization_monitor
    self.model = model
    self.log = log
    print("Minimizing using reference map...", file=self.log)
    self.log.flush()

    # copy-paste from cctbx_project/mmtbx/refinement/geometry_minimization.py:
    # minimize_wrapper_for_ramachandran
    self.model.get_restraints_manager().geometry.pair_proxies(
        sites_cart=self.model.get_sites_cart())

    ncs_restraints_group_list = self.model.get_ncs_groups()
    if ncs_restraints_group_list is None:
      ncs_restraints_group_list = []
    ncs_groups=None
    if len(ncs_restraints_group_list) > 0:
      ncs_groups=ncs_restraints_group_list

    min_monitor = minimization_monitor(
        number_of_cycles=number_of_cycles,
        max_number_of_cycles=20,
        cycles_to_converge=cycles_to_converge,
        mode=min_mode)
    selection_real_space = None
    import mmtbx.refinement.real_space.weight
    self.w = 1
    print("number_of_cycles", number_of_cycles, file=log)
    print("Stats before minimization:", file=log)
    ms = self.model.geometry_statistics()
    ms.show(log=log)

    while min_monitor.need_more_cycles():
      print("Cycle number", min_monitor.get_current_cycle_n(), file=self.log)
      self.model.get_restraints_manager().geometry.\
          update_ramachandran_restraints_phi_psi_targets(
              hierarchy=self.model.get_hierarchy())
      print("  Updating rotamer restraints...", file=self.log)
      add_rotamer_restraints(
        pdb_hierarchy      = self.model.get_hierarchy(),
        restraints_manager = self.model.get_restraints_manager(),
        selection          = None,
        sigma              = 5,
        mode               = "fix_outliers",
        accept_allowed     = False,
        mon_lib_srv        = self.model.get_mon_lib_srv(),
        rotamer_manager    = self.model.get_rotamer_manager())
      self.model.set_sites_cart_from_hierarchy()

      if min_monitor.need_weight_optimization():
        # if self.w is None:
        print("  Determining weight...", file=self.log)
        self.log.flush()
        self.weight = mmtbx.refinement.real_space.weight.run(
            map_data                    = target_map,
            xray_structure              = self.model.get_xray_structure(),
            pdb_hierarchy               = self.model.get_hierarchy(),
            geometry_restraints_manager = self.model.get_restraints_manager(),
            rms_bonds_limit             = 0.015,
            rms_angles_limit            = 1.0,
            ncs_groups                  = ncs_groups)

        # division is to put more weight onto restraints. Checked. Works.
        self.w = self.weight.weight/3.0
        # self.w = self.weight.weight/15.0
        # self.w = 0
        # self.w = self.weight.weight
        for s in self.weight.msg_strings:
          print(s, file=self.log)

      if ncs_restraints_group_list is None or len(ncs_restraints_group_list)==0:
        #No NCS
        print("  Minimizing...", file=self.log)
        print("     with weight %f" % self.w, file=self.log)
        self.log.flush()
        refine_object = simple(
            target_map                  = target_map,
            selection                   = None,
            max_iterations              = 150,
            geometry_restraints_manager = self.model.get_restraints_manager().geometry,
            selection_real_space        = selection_real_space,
            states_accumulator          = None,
            ncs_groups                  = ncs_groups)
        refine_object.refine(weight = self.w, xray_structure = self.model.get_xray_structure())
        self.rmsd_bonds_final, self.rmsd_angles_final = refine_object.rmsds()
        print("RMSDS:", self.rmsd_bonds_final, self.rmsd_angles_final, file=log)
        # print >> log, "sizes:", len(refine_object.sites_cart()), len(self.xrs.scatterers())
        self.model.set_sites_cart(refine_object.sites_cart(), update_grm=True)
        # print >> log, "sizes", self.xrs.scatterers()
      else:
        # Yes NCS
        # copy-paste from macro_cycle_real_space.py
        # !!! Don't rearrange NCS groups here because master was just fixed!
        # import mmtbx.ncs.ncs_utils as nu
        # nu.get_list_of_best_ncs_copy_map_correlation(
        #     ncs_groups     = ncs_restraints_group_list,
        #     xray_structure = self.model.get_xray_structure(),
        #     map_data       = target_map,
        #     d_min          = 3)
        print("  Minimizing... (NCS)", file=self.log)
        tfg_obj = mmtbx.refinement.minimization_ncs_constraints.\
          target_function_and_grads_real_space(
            map_data                   = target_map,
            xray_structure             = self.model.get_xray_structure(),
            ncs_restraints_group_list  = ncs_restraints_group_list,
            refine_selection           = None,
            real_space_gradients_delta = 1,
            restraints_manager         = self.model.get_restraints_manager(),
            data_weight                = self.w,
            refine_sites               = True)
        minimized = mmtbx.refinement.minimization_ncs_constraints.lbfgs(
          target_and_grads_object      = tfg_obj,
          xray_structure               = self.model.get_xray_structure(),
          ncs_restraints_group_list    = ncs_restraints_group_list,
          refine_selection             = None,
          finite_grad_differences_test = False,
          max_iterations               = 100,
          refine_sites                 = True)
        self.model.set_sites_cart(tfg_obj.xray_structure.sites_cart())
      ncs_restraints_group_list.recalculate_ncs_transforms(self.model.get_sites_cart())
      ms = self.model.geometry_statistics()
      min_monitor.save_cycle_results(geometry=ms)
      ms.show(log=log)
