from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.lbfgs
from scitbx import matrix
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args
import sys
import boost_adaptbx.boost.python as bp
from six.moves import range
cctbx_maptbx_ext = bp.import_ext("cctbx_maptbx_ext")
from cctbx import miller
import mmtbx.utils

def real_space_rigid_body_gradients_simple(
      unit_cell,
      density_map,
      sites_cart_0,
      center_of_mass,
      q,
      unit_quaternion_delta=0.01,
      translation_delta=0.3):
  result = flex.double()
  q_delta = q.deep_copy()
  def get(i, delta):
    fs = []
    for signed_delta in [delta, -delta]:
      q_delta[i] = q[i] + signed_delta
      aja = matrix.rt(scitbx.rigid_body.joint_lib_six_dof_aja_simplified(
        center_of_mass=center_of_mass,
        q=q_delta))
      sites_cart_delta = aja * sites_cart_0
      rs_f = maptbx.real_space_target_simple(
        unit_cell=unit_cell,
        density_map=density_map,
        sites_cart=sites_cart_delta,
        selection=flex.bool(sites_cart_delta.size(),True))
      fs.append(rs_f)
    result.append((fs[0]-fs[1])/(2*delta))
  for i in range(4): get(i=i, delta=unit_quaternion_delta)
  for i in range(3): get(i=i+4, delta=translation_delta)
  return result

class refine2(object):

  def __init__(self,
        sites_cart,
        density_map,
        real_space_gradients_delta,
        lbfgs_termination_params,
        unit_cell,
        states_collector=None):
    self.states_collector = states_collector
    self.density_map = density_map
    self.real_space_gradients_delta = real_space_gradients_delta
    #
    self.unit_cell = unit_cell
    self.sites_cart = sites_cart
    self.center_of_mass = self.sites_cart.mean()
    tardy_tree = scitbx.graph.tardy_tree.construct(
      n_vertices=self.sites_cart.size(),
      edge_list="all_in_one_rigid_body") \
        .build_tree() \
        .fix_near_singular_hinges(sites=None)
    self.tardy_model = scitbx.rigid_body.tardy_model(
      labels=None,
      sites=self.sites_cart,
      masses=flex.double(self.sites_cart.size(), 1),
      tardy_tree=tardy_tree,
      potential_obj=self)
    self.x = self.tardy_model.pack_q()
    assert self.x.size() == 7 # other cases not implemented
    #
    self.number_of_function_evaluations = -1
    self.f_start, self.g_start = self.compute_functional_and_gradients()
    self.rs_f_start = self.rs_f
    lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      exception_handling_params = lbfgs_exception_handling_params)
    self.f_final, self.g_final = self.compute_functional_and_gradients()
    self.rs_f_final = self.rs_f
    del self.rs_f
    del self.x
    del self.center_of_mass
    del self.unit_cell

  def compute_functional_and_gradients(self):
    if (self.number_of_function_evaluations == 0):
      self.number_of_function_evaluations += 1
      return self.f_start, self.g_start
    self.number_of_function_evaluations += 1
    self.tardy_model.unpack_q(q_packed=self.x)
    self.sites_cart = self.tardy_model.sites_moved()
    if(self.states_collector is not None):
      self.states_collector.add(sites_cart = self.sites_cart)
    rs_f = maptbx.real_space_target_simple(
      unit_cell=self.unit_cell,
      density_map=self.density_map,
      sites_cart=self.sites_cart,
      selection=flex.bool(self.sites_cart.size(),True))
    rs_g = real_space_rigid_body_gradients_simple(
      unit_cell=self.unit_cell,
      density_map=self.density_map,
      sites_cart_0=self.sites_cart,
      center_of_mass=self.center_of_mass,
      q=self.x)
    self.rs_f = rs_f
    f = -1*rs_f
    g = -1*rs_g
    return f, g.as_double()

  def d_e_pot_d_sites(self, sites_moved):
    result = self.__d_e_pot_d_sites
    del self.__d_e_pot_d_sites
    return result

class refine(object):

  def __init__(self,
        residue,
        density_map,
        geometry_restraints_manager,
        real_space_target_weight,
        real_space_gradients_delta,
        lbfgs_termination_params,
        unit_cell,
        cctbx_geometry_restraints_flags=None,
        states_collector=None):
    self.states_collector = states_collector
    self.cctbx_geometry_restraints_flags = cctbx_geometry_restraints_flags
    self.residue = residue
    self.density_map = density_map
    self.geometry_restraints_manager = geometry_restraints_manager
    self.real_space_gradients_delta = real_space_gradients_delta
    self.real_space_target_weight = real_space_target_weight
    #
    self.unit_cell = unit_cell
    self.sites_cart_residue_0 = residue.atoms().extract_xyz()
    self.residue_center_of_mass = self.sites_cart_residue_0.mean()
    residue_tardy_tree = scitbx.graph.tardy_tree.construct(
      n_vertices=self.sites_cart_residue_0.size(),
      edge_list="all_in_one_rigid_body") \
        .build_tree() \
        .fix_near_singular_hinges(sites=None)
    self.residue_tardy_model = scitbx.rigid_body.tardy_model(
      labels=None,
      sites=self.sites_cart_residue_0,
      masses=flex.double(self.sites_cart_residue_0.size(), 1),
      tardy_tree=residue_tardy_tree,
      potential_obj=self)
    self.x = self.residue_tardy_model.pack_q()
    assert self.x.size() == 7 # other cases not implemented
    #
    self.number_of_function_evaluations = -1
    self.f_start, self.g_start = self.compute_functional_and_gradients()
    self.rs_f_start = self.rs_f
    lbfgs_exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
      ignore_line_search_failed_step_at_lower_bound = True,
      ignore_line_search_failed_step_at_upper_bound = True,
      ignore_line_search_failed_maxfev              = True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs_termination_params,
      exception_handling_params = lbfgs_exception_handling_params)
    self.f_final, self.g_final = self.compute_functional_and_gradients()
    self.rs_f_final = self.rs_f
    del self.rs_f
    del self.x
    del self.residue_center_of_mass
    del self.sites_cart_residue_0
    del self.unit_cell

  def compute_functional_and_gradients(self):
    if (self.number_of_function_evaluations == 0):
      self.number_of_function_evaluations += 1
      return self.f_start, self.g_start
    self.number_of_function_evaluations += 1
    self.residue_tardy_model.unpack_q(q_packed=self.x)
    self.sites_cart_residue = self.residue_tardy_model.sites_moved()
    if(self.states_collector is not None):
      self.states_collector.add(sites_cart = self.sites_cart_residue)
    rs_f = maptbx.real_space_target_simple(
      unit_cell=self.unit_cell,
      density_map=self.density_map,
      sites_cart=self.sites_cart_residue,
      selection=flex.bool(self.sites_cart_residue.size(),True))
    rs_g = real_space_rigid_body_gradients_simple(
      unit_cell=self.unit_cell,
      density_map=self.density_map,
      sites_cart_0=self.sites_cart_residue_0,
      center_of_mass=self.residue_center_of_mass,
      q=self.x)
    self.rs_f = rs_f
    rs_f *= -self.real_space_target_weight
    rs_g *= -self.real_space_target_weight
    if (self.geometry_restraints_manager is None):
      f = rs_f
      g = rs_g
    else:
      gr_e = self.geometry_restraints_manager.energies_sites(
        sites_cart=self.sites_cart_residue,
        flags = self.cctbx_geometry_restraints_flags,
        compute_gradients=True)
      self.__d_e_pot_d_sites = gr_e.gradients
      f = rs_f + gr_e.target
      g = rs_g + self.residue_tardy_model.d_e_pot_d_q_packed()
    return f, g.as_double()

  def d_e_pot_d_sites(self, sites_moved):
    result = self.__d_e_pot_d_sites
    del self.__d_e_pot_d_sites
    return result

class refine_mz(object):
  """
  Efficient real-space rigid-body refinement. Analog of MZ rigid-body refinement
  in reciprocal space. Whole content of pdb_hierarchy is treated as one rigid
  group.
  """

  def __init__(
        self,
        map_data,
        pdb_hierarchy,  # XXX redundant inputs
        xray_structure, # XXX redundant inputs
        d_min,
        use_mask=False,
        masking_atom_radius=5,
        max_iterations=50,
        macro_cycles=1,
        use_one_zone=False,
        prefix="",
        log=None):
    adopt_init_args(self, locals())
    self.cc_best = None
    self.sites_cart_best = None
    if(self.log is None): self.log = sys.stdout
    self.sites_cart_start = self.xray_structure.sites_cart()
    assert approx_equal(self.pdb_hierarchy.atoms().extract_xyz(),
      self.sites_cart_start, 1.e-3)
    self.crystal_gridding = maptbx.crystal_gridding(
      unit_cell             = self.xray_structure.unit_cell(),
      space_group_info      = self.xray_structure.space_group_info(),
      pre_determined_n_real = self.map_data.all())
    self.complete_set = miller.build_set(
      crystal_symmetry = self.xray_structure.crystal_symmetry(),
      anomalous_flag   = False,
      d_min            = self.d_min)
    self.cc_start = self._get_cc()
    self._show_and_track(d_min = self.d_min)
    if not use_one_zone:
      self.d_mins = self._get_mz_resolution_limits()
    else:
      self.d_mins = [d_min, ]
    for macro_cycle in range(self.macro_cycles):
      self._refine()
    self.xray_structure.set_sites_cart(self.sites_cart_best)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

  def _show_and_track(self, d_min):
    cc = self._get_cc()
    s2 = self.xray_structure.sites_cart()
    if(self.cc_best is None or cc>self.cc_best):
      self.cc_best = cc
      self.sites_cart_best = s2.deep_copy()
    if(self.log):
      fmt="%sd_min=%6.2f CC=%6.4f (best to keep CC=%6.4f), moved from start (max/mean)=%s"
      s1 = self.sites_cart_start
      d = flex.sqrt((s1-s2).dot()).min_max_mean().as_tuple()[1:]
      d_str = "%6.3f %6.3f"%d
      print(fmt%(self.prefix, d_min, cc, self.cc_best, d_str), file=self.log)
      revert = (cc<self.cc_start and abs(cc-self.cc_start)>0.1) or \
               d[0] > 10 and (cc<self.cc_start or abs(cc-self.cc_start)<0.1) or \
               cc < 0.1
      if(revert):
        self.sites_cart_best = self.sites_cart_start
        self.xray_structure.set_sites_cart(self.sites_cart_start)
        self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
        print("   >>> reversed", d[0], cc, self.cc_start, file=self.log)

  def _refine(self):
    self.lbfgs_termination_params=scitbx.lbfgs.termination_parameters(
      max_iterations = self.max_iterations)
    mask_data = self._get_mask()
    for d_min in self.d_mins:
      md = self._get_map_at_d_min(d_min=d_min)
      if(mask_data is not None):
        md = md*mask_data.as_double()
      minimized = refine(
        residue                     = self.pdb_hierarchy,
        density_map                 = md,
        geometry_restraints_manager = None,
        real_space_target_weight    = 1,
        real_space_gradients_delta  = d_min*0.25,
        lbfgs_termination_params    = self.lbfgs_termination_params,
        unit_cell                   = self.xray_structure.unit_cell())
      self.xray_structure.set_sites_cart(minimized.sites_cart_residue)
      self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)
      self._show_and_track(d_min = d_min)

  def _get_mz_resolution_limits(self):
    # lowest resolution: first zone
    n_ref_lowest = 0
    d_spacing = self.complete_set.d_spacings().sort().data()
    if(d_spacing.size()<500): return ( self.complete_set.d_min(), )
    d_0 = min(d_spacing[499], 15)
    if(d_0>8.0): d_1 = 8.0
    else:
      return d_0, self.d_min
    # second zone
    if(d_1>4.0 and self.d_min<4.0):
      d_2 = 4.0
    else:
      d_2 = max(self.d_min, 4.0)
    return d_0, d_1, d_2

  def _get_mask(self):
    mask_data = None
    if(self.use_mask):
      xrs_p1 = self.xray_structure.expand_to_p1(sites_mod_positive=True)
      radii = flex.double(xrs_p1.scatterers().size(), self.masking_atom_radius)
      mask_data = cctbx_maptbx_ext.mask(
        sites_frac                  = xrs_p1.sites_frac(),
        unit_cell                   = xrs_p1.unit_cell(),
        n_real                      = self.map_data.all(),
        mask_value_inside_molecule  = 1,
        mask_value_outside_molecule = 0,
        radii                       = radii)
    return mask_data

  def _get_map(self):
    f_calc = self.xray_structure.structure_factors(d_min=self.d_min).f_calc()
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = f_calc)
    fft_map.apply_sigma_scaling() # XXX not really needed
    return fft_map.real_map_unpadded()

  def _get_cc(self):
    return flex.linear_correlation(
      x=self.map_data.as_1d(),
      y=self._get_map().as_1d()).coefficient()

  def _get_map_at_d_min(self, d_min):
    done = False
    cntr = 0
    while not done:
      if(cntr>50):
        raise RuntimeError("Number of trial resolution increments exceeded.")
      try:
        f_obs_cmpl = self.complete_set.resolution_filter(
          d_min=d_min).structure_factors_from_map(
            map            = self.map_data,
            use_scale      = True,
            anomalous_flag = False,
            use_sg         = True)
        done = True
      except KeyboardInterrupt: raise
      except Exception as e:
        if(str(e)=="cctbx Error: Miller index not in structure factor map."):
          d_min += 0.25
      cntr+=1
    fft_map = miller.fft_map(
      crystal_gridding     = self.crystal_gridding,
      fourier_coefficients = f_obs_cmpl)
    fft_map.apply_sigma_scaling()
    return fft_map.real_map_unpadded()

class refine_groups(object):

  def __init__(
        self,
        map_data,
        pdb_hierarchy,
        xray_structure,
        macro_cycles,
        d_min):
    self.pdb_hierarchy = pdb_hierarchy
    self.xray_structure = xray_structure
    # sanity check
    sites_cart_result = self.xray_structure.sites_cart().deep_copy()
    assert approx_equal(sites_cart_result,
      self.pdb_hierarchy.atoms().extract_xyz(), 1.e-4)
    #
    cs = self.xray_structure.crystal_symmetry()
    for i_chain, chain in enumerate(self.pdb_hierarchy.chains()):
      print("chain:", chain.id)
      selection = chain.atoms().extract_i_seq()
      ph = pdb_hierarchy.select(selection)
      #
      xrs_tmp = xray_structure.select(selection)
      box = mmtbx.utils.extract_box_around_model_and_map(
        xray_structure = xrs_tmp,
        map_data       = map_data,
        box_cushion    = 5,
       )
      #
      shift_back = [-box.shift_cart[i] for i in range(3)]
      ph_b       = box.pdb_hierarchy_box
      md_b       = box.map_box
      xrs_b      = box.xray_structure_box
      #
      minimized = refine_mz(
        map_data       = md_b,
        pdb_hierarchy  = ph_b,
        xray_structure = xrs_b,
        d_min          = d_min,
        macro_cycles   = macro_cycles,
        log            = None,
        prefix="  ")
      sites_cart_result = sites_cart_result.set_selected(
        selection, minimized.sites_cart_best+shift_back)
    self.xray_structure.set_sites_cart(sites_cart_result)
    self.pdb_hierarchy.adopt_xray_structure(self.xray_structure)

def translation_search(model, map_data, shifts=None):
  """
  Rigid-body translation grid search in real-space
  """
  sites_cart_best = model.get_sites_cart()
  cs = model.crystal_symmetry()
  fm = cs.unit_cell().fractionalization_matrix()
  if shifts is None:
    shifts = [i/10. for i in range(0,101, 2)]
  t_best = -999999999.
  sites_cart = sites_cart_best.deep_copy()
  size = sites_cart.size()
  best_shift = None
  for sh in shifts:
    for x in [-1,0,1]:
      for y in [-1,0,1]:
        for z in [-1,0,1]:
          sh_ = [x*sh,y*sh,z*sh]
          sites_cart_shifted = sites_cart + \
            flex.vec3_double(sites_cart.size(), sh_)
          sites_frac_shifted = fm*sites_cart_shifted
          t = maptbx.map_sum_at_sites_frac(
            map_data   = map_data,
            sites_frac = sites_frac_shifted)/size
          if(t>t_best):
            t_best=t
            best_shift = sh_[:]
            sites_cart_best=sites_cart_shifted.deep_copy()
  #print("best_shift:", best_shift)
  model.set_sites_cart(sites_cart = sites_cart_best)
  return best_shift

def protocol_1(model, map_data, macro_cycles=6, max_iterations=25):
  """
  Rigid-body grid search translation search combined with rigid body gradient
  based minimization in real-space
  """
  shift = translation_search(model=model, map_data=map_data)
  #
  for it in range(macro_cycles):
    delta = 0.01
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    o = refine2(
      sites_cart                      = model.get_sites_cart(),
      density_map                     = map_data,
      real_space_gradients_delta      = delta,
      lbfgs_termination_params        = lbfgs_termination_params,
      unit_cell                       = model.crystal_symmetry().unit_cell(),
      states_collector                = None)
    model.set_sites_cart(sites_cart = o.sites_cart)
    shift = translation_search(model=model, map_data=map_data)
    max_shift = flex.max(flex.abs(flex.double(shift)))
    if max_shift < 0.09: break
  #
  shifts = [x * 0.01 for x in range(0, 21)]
  shift = translation_search(model=model, map_data=map_data, shifts=shifts)
