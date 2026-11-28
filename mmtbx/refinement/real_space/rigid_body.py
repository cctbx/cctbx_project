from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx.array_family import flex
import scitbx.rigid_body
import scitbx.graph.tardy_tree
import scitbx.lbfgs
from scitbx import matrix
from scitbx.math import superpose
from libtbx.test_utils import approx_equal
from libtbx import adopt_init_args
import sys, math
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
    return mz_resolution_limits(
      complete_set = self.complete_set, d_min = self.d_min)

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
    return map_at_d_min(
      map_data         = self.map_data,
      complete_set     = self.complete_set,
      crystal_gridding = self.crystal_gridding,
      d_min            = d_min)

def mz_resolution_limits(complete_set, d_min):
  """
  Resolution limits of zones for multi-zone rigid-body refinement, from the
  lowest resolution (largest radius of convergence) to the highest.
  """
  d_spacing = complete_set.d_spacings().sort().data()
  if(d_spacing.size()<500): return ( complete_set.d_min(), )
  d_0 = min(d_spacing[499], 15)
  if(d_0>8.0): d_1 = 8.0
  else:
    return d_0, d_min
  # second zone
  if(d_1>4.0 and d_min<4.0):
    d_2 = 4.0
  else:
    d_2 = max(d_min, 4.0)
  return d_0, d_1, d_2

def map_at_d_min(map_data, complete_set, crystal_gridding, d_min):
  """
  Map with Fourier coefficients beyond d_min removed (sigma-scaled).
  """
  done = False
  cntr = 0
  while not done:
    if(cntr>50):
      raise RuntimeError("Number of trial resolution increments exceeded.")
    try:
      f_obs_cmpl = complete_set.resolution_filter(
        d_min=d_min).structure_factors_from_map(
          map            = map_data,
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
    crystal_gridding     = crystal_gridding,
    fourier_coefficients = f_obs_cmpl)
  fft_map.apply_sigma_scaling()
  return fft_map.real_map_unpadded()

def _subset(size, max_sites):
  """
  Selection of every n-th of size atoms so that about max_sites are selected;
  None if all are to be used.
  """
  if(max_sites is None or size <= max_sites): return None
  return flex.size_t(range(0, size, int(math.ceil(size/max_sites))))

def translation_search(model, map_data, shifts=None, max_shift=10.,
                       max_sites=None):
  """
  Rigid-body translation grid search in real-space: the model is moved along
  the 26 directions of a 3x3x3 grid by each of the shifts (default: 0 to
  max_shift, step 0.2 A) and left at the position with the largest mean map
  value at atom centers. Returns the shift applied.
  max_sites: if given and the model has more atoms, every n-th atom is used
  so that about max_sites atoms are used (a few thousand atoms are enough to
  rank shifts and it is much faster).
  """
  sites_cart = model.get_sites_cart()
  fm = matrix.sqr(
    model.crystal_symmetry().unit_cell().fractionalization_matrix())
  sites_frac = model.get_sites_frac()
  subset = _subset(sites_frac.size(), max_sites)
  if(subset is not None): sites_frac = sites_frac.select(subset)
  size = sites_frac.size()
  if shifts is None:
    shifts = [i/10. for i in range(0, int(round(max_shift*10))+1, 2)]
  t_best = -999999999.
  best_shift = [0,0,0]
  tried = set()
  for sh in shifts:
    for x in [-1,0,1]:
      for y in [-1,0,1]:
        for z in [-1,0,1]:
          sh_ = [x*sh,y*sh,z*sh]
          if(tuple(sh_) in tried): continue # zero shift comes up many times
          tried.add(tuple(sh_))
          t = maptbx.map_sum_at_sites_frac(
            map_data   = map_data,
            sites_frac = sites_frac + (fm*matrix.col(sh_)).elems)/size
          if(t>t_best):
            t_best=t
            best_shift = sh_[:]
  #print("best_shift:", best_shift)
  model.set_sites_cart(sites_cart = sites_cart + tuple(best_shift))
  return best_shift

def protocol_1(model, map_data, macro_cycles=6, max_iterations=25,
               max_shift=10., max_sites=None):
  """
  Rigid-body grid search translation search combined with rigid body gradient
  based minimization in real-space. The whole model is one rigid body; the
  target is the mean map value at atom centers. max_shift: range of the
  translation grid search (A). max_sites: if given and the model has more
  atoms, position and orientation are determined using every n-th atom (see
  translation_search) and applied to all atoms.
  Used by phenix.magref (whole model) and, through refine_groups below, by
  phenix.real_space_refine.
  """
  shift = translation_search(model=model, map_data=map_data,
    max_shift=max_shift, max_sites=max_sites)
  if(model.size() < 3): return # too few atoms to define orientation
  #
  for it in range(macro_cycles):
    delta = 0.01
    lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
      max_iterations = max_iterations)
    sites_cart = model.get_sites_cart()
    subset = _subset(sites_cart.size(), max_sites)
    sites_cart_used = sites_cart
    if(subset is not None): sites_cart_used = sites_cart.select(subset)
    o = refine2(
      sites_cart                      = sites_cart_used,
      density_map                     = map_data,
      real_space_gradients_delta      = delta,
      lbfgs_termination_params        = lbfgs_termination_params,
      unit_cell                       = model.crystal_symmetry().unit_cell(),
      states_collector                = None)
    # Keep the result of minimization only if the target (sum of map values
    # at atom centers) did not get worse: a failed line search can send the
    # model far away.
    if(o.rs_f_final >= o.rs_f_start):
      if(subset is None):
        model.set_sites_cart(sites_cart = o.sites_cart)
      else: # apply the rigid-body move of the atoms used to all atoms
        lsq = superpose.least_squares_fit(
          reference_sites = o.sites_cart,
          other_sites     = sites_cart_used)
        model.set_sites_cart(
          sites_cart = lsq.r.elems * sites_cart + lsq.t.elems)
    shift = translation_search(model=model, map_data=map_data,
      max_shift=max_shift, max_sites=max_sites)
    max_shift_found = flex.max(flex.abs(flex.double(shift)))
    if max_shift_found < 0.09: break
  #
  shifts = [x * 0.01 for x in range(0, 21)]
  shift = translation_search(model=model, map_data=map_data, shifts=shifts,
    max_sites=max_sites)

def mean_map_value(model, map_data):
  """
  Mean map value at atom centers: the target of protocol_1.
  """
  return maptbx.map_sum_at_sites_frac(
    map_data   = map_data,
    sites_frac = model.get_sites_frac())/model.size()

def _refine_in_box(args):
  """
  Refine one group in its box with each of the search ranges. Returns
  group index and a list of (mean map value, sites_cart), the first one for
  the start position.
  """
  i, box_mmm, local_shift, max_shift, d_min, max_sites = args
  model = box_mmm.model()
  map_data = box_mmm.map_manager().map_data()
  sites_cart_start = model.get_sites_cart().deep_copy()
  result = [(mean_map_value(model, map_data), sites_cart_start)]
  def add():
    result.append(
      (mean_map_value(model, map_data), model.get_sites_cart().deep_copy()))
  # local
  protocol_1(model = model, map_data = map_data, max_shift = local_shift,
    max_sites = max_sites)
  add()
  # long-range
  if(max_shift is not None and max_shift > local_shift):
    model.set_sites_cart(sites_cart_start)
    if(d_min is not None):
      # multi-zone: maps of increasing resolution; the lowest resolution
      # gives the largest radius of convergence
      cs = model.crystal_symmetry()
      crystal_gridding = maptbx.crystal_gridding(
        unit_cell             = cs.unit_cell(),
        space_group_info      = cs.space_group_info(),
        pre_determined_n_real = map_data.all())
      complete_set = miller.build_set(
        crystal_symmetry = cs,
        anomalous_flag   = False,
        d_min            = d_min)
      for d in mz_resolution_limits(complete_set = complete_set, d_min = d_min):
        protocol_1(
          model     = model,
          map_data  = map_at_d_min(
            map_data         = map_data,
            complete_set     = complete_set,
            crystal_gridding = crystal_gridding,
            d_min            = d),
          max_shift = max_shift,
          max_sites = max_sites)
      protocol_1(model = model, map_data = map_data, max_shift = local_shift,
        max_sites = max_sites)
    else:
      protocol_1(model = model, map_data = map_data, max_shift = max_shift,
        max_sites = max_sites)
    add()
  model.set_sites_cart(sites_cart_start)
  return i, result

class refine_groups(object):
  """
  Real-space rigid-body refinement of several rigid groups (selections) of
  the model in map_model_manager, each with protocol_1 in a box cut out of
  the map around the group. Model in map_model_manager is updated.

  Groups are refined independently of each other, and density is not removed
  or masked around the other groups (masking biases groups that are in
  contact and leaves nothing to fit in tightly packed assemblies). The price
  is that a group could move onto density that is not its own: of a
  neighbour, or of a part that is not modelled. In assemblies of identical
  molecules such a place is as good as its own (for example the next layer
  of a fibril, 4.8 A away). Therefore each group is refined twice from where
  it is:
    - local: short-range translation search (local_shift, 1 A). A group
      displaced by more than that still gets to its place in several steps
      (grid search and minimization alternate) as long as the target
      improves along the way; expect convergence from about 2 A / 5 degrees;
    - long-range: translation search up to max_shift (5 A); if d_min (map
      resolution) is given, first in maps truncated to lower resolution
      (multi-zone, the lowest resolution gives the largest radius of
      convergence: tens of A for large groups), then local as above;
  and the long-range result is taken only if its mean map value at atom
  centers is higher than that of the local result by more than the fraction
  min_gain (a long move must pay; a move to an equivalent place does not).
  The local result is taken only if it is better than the start position by
  more than the fraction min_gain_local: a group that is in place is not
  moved (the optimum of the target is within about 0.1-0.3 A of the best
  position, depending on resolution, and such a shift gains less than that).
  A selection that selects nothing is skipped. Groups of fewer than
  min_size_long_range atoms are refined locally only (for a few atoms any
  strong peak nearby "pays").
  A result is not used if the group moved (mean over atoms) by more than
  max_move (if given), or if it makes the group overlap with atoms outside
  the group:
  the number of such atoms within overlap_distance of atoms of the group,
  relative to group size, has grown by more than overlap_tolerance (the
  distance is short, 1 A, so that contacts, hydrogen bonds and a few covalent
  links between groups do not count: this is about one group sitting on top
  of another). Then the next best result is used; the last one is the start
  position.

  nproc: groups are done in parallel; the result does not depend on nproc.
  """

  def __init__(
        self,
        map_model_manager,
        selections,
        d_min=None,
        local_shift=1.,
        max_shift=5.,
        min_gain=0.1,
        min_gain_local=0.01,
        max_move=None,
        overlap_distance=1.0,
        overlap_tolerance=0.05,
        max_sites=3000,
        min_size_long_range=10,
        nproc=1,
        log=None):
    self.mmm = map_model_manager
    self.selections = selections
    self.overlap_distance = overlap_distance
    self.log = log
    model = self.mmm.model()
    size = model.size()
    for sel in selections:
      assert isinstance(sel, flex.bool) and sel.size() == size
    sites_cart_start = model.get_sites_cart().deep_copy()
    overlap_start = [self._overlap(model, sel) for sel in selections]
    # empty selections get no box: keep track of which group each box is
    not_empty = [i for i, sel in enumerate(selections) if sel.count(True) > 0]
    argss = []
    if(len(not_empty) > 0):
      box_info = self.mmm.split_up_map_and_model_by_supplied_selections(
        selection_list = [selections[i] for i in not_empty])
      assert len(box_info.mmm_list) == len(not_empty)
      for i, box_mmm in zip(not_empty, box_info.mmm_list):
        ms = max_shift
        if(selections[i].count(True) < min_size_long_range): ms = None
        argss.append([i, box_mmm, local_shift, ms, d_min, max_sites])
    if(nproc > 1 and len(argss) > 1):
      from libtbx import easy_mp
      results = easy_mp.pool_map(
        processes  = nproc,
        fixed_func = _refine_in_box,
        args       = argss)
    else:
      results = [_refine_in_box(args) for args in argss]
    # Candidate positions of each group in order of preference, as shifts of
    # atoms from the start position (boxes are shifted copies of the model).
    self.candidates = [
      [["not moved", 0, flex.vec3_double(sel.count(True), (0,0,0))]]
      for sel in selections]
    for i, result in results:
      score_start, sites_start = result[0]
      trials = [] # (label, score, shifts)
      for label, (score, sites) in zip(["local", "long-range"], result[1:]):
        d = sites - sites_start
        if(max_move is not None and
           flex.mean(flex.sqrt(d.dot())) > max_move): continue
        trials.append([label, score, d])
      if(len(trials) == 2 and
         trials[1][1] - trials[0][1] > min_gain*abs(trials[0][1])):
        trials.reverse()
      elif(len(trials) == 2):
        trials = trials[:1]
      trials = [t for t in trials
        if t[1] - score_start > min_gain_local*abs(score_start)]
      trials.append(["not moved", score_start,
        flex.vec3_double(sites_start.size(), (0,0,0))])
      self.candidates[i] = trials
    choice = [0,]*len(selections)
    def set_sites():
      sites_cart = sites_cart_start.deep_copy()
      for i, sel in enumerate(selections):
        sites_cart.set_selected(sel,
          sites_cart_start.select(sel) + self.candidates[i][choice[i]][2])
      model.set_sites_cart(sites_cart)
      return sites_cart
    # Give up results that end in overlap with other atoms: the group that
    # moved most goes first (it is the one that came to the others).
    while True:
      sites_cart = set_sites()
      moved = flex.sqrt((sites_cart - sites_cart_start).dot())
      overlapping = []
      for i, sel in enumerate(selections):
        if(choice[i] == len(self.candidates[i])-1): continue # at start
        if(self._overlap(model, sel) > overlap_start[i] + overlap_tolerance):
          overlapping.append([flex.mean(moved.select(sel)), i])
      if(len(overlapping) == 0): break
      overlapping.sort()
      d, i = overlapping[-1]
      self._print("  group %d: %s result (moved %.2f A) runs into other atoms,"
        " not used"%(i, self.candidates[i][choice[i]][0], d))
      choice[i] += 1
    sites_cart = model.get_sites_cart()
    moved = flex.sqrt((sites_cart - sites_cart_start).dot())
    self.moved = [flex.mean_default(moved.select(sel), 0) for sel in selections]
    self.result = [self.candidates[i][choice[i]][0]
                   for i in range(len(selections))]
    for i, d in enumerate(self.moved):
      c = self.candidates[i]
      self._print(
        "  group %d: %-10s moved (mean): %6.3f A  mean map value: %s -> %s"%(
        i, self.result[i], d, "%.4g"%c[-1][1], "%.4g"%c[choice[i]][1]))

  def _overlap(self, model, selection):
    """
    Number of atoms outside the group within overlap_distance of atoms of
    the group, divided by the number of atoms in the group.
    """
    n = selection.count(True)
    if(n == 0 or n == selection.size()): return 0
    near = model.get_xray_structure().selection_within(
      radius    = self.overlap_distance,
      selection = selection)
    return (near & ~selection).count(True)/n

  def _print(self, m):
    if(self.log is not None): print(m, file=self.log)
