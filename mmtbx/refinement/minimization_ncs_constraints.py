from __future__ import division
import cctbx.maptbx.real_space_refinement_simple
from scitbx.array_family import flex
from libtbx import adopt_init_args
import mmtbx.utils.ncs_utils as nu
from cctbx import xray
import scitbx.lbfgs

def grads_asu_to_one_ncs(
      ncs_restraints_group_list,
      extended_ncs_selection,
      grad,
      refine_sites):
  """
  Apply NCS constraints to reduce the gradients corresponding to whole ASU, to
  gradients corresponding to one NCS copy.

  Arguments:
  ncs_restraints_group_list: list of ncs_restraint_group objects
  grad: gradient of the complete asu
  refine_sites: (bool) Flag indicating sites refinement
  """
  # TODO: write a test for this function
  for nrg in ncs_restraints_group_list:
    ncs_selection = nrg.master_iselection
    for ncs_copy in nrg.copies:
      asu_selection = ncs_copy.copy_iselection
      g = grad.select(asu_selection)
      if(refine_sites):
        # apply inverse transformation
        rt = ncs_copy.r.transpose().elems
        g = (rt*g)
      grad.set_selected(ncs_selection, g + grad.select(ncs_selection))

  g_ncs = grad.select(extended_ncs_selection)
  if(refine_sites): g_ncs = flex.vec3_double(g_ncs)
  assert type(grad)==type(g_ncs)
  return g_ncs

def grads_one_ncs_to_asu(ncs_restraints_group_list,
                         total_asu_length,
                         extended_ncs_selection,
                         ncs_grad):
  """
  Expand average gradient of a single NCS to all ASU
  (only for u_iso refinement)

  Arguments:
  ncs_restraints_group_list: list of ncs_restraint_group objects
  ncs_grad: gradient of ta single ncs
  """
  # TODO: Add test to make sure this works, check asu_to_ncs_map
  g_length = total_asu_length
  if isinstance(ncs_grad,flex.vec3_double):
    g = flex.vec3_double([(0.0,0.0,0.0)]*g_length)
  elif isinstance(ncs_grad,flex.double):
    g = flex.double([0.0]*g_length)
  else:
    raise TypeError('Non supported grad type')
  # update newly created flex.vec3 with master NCS info
  g.set_selected(extended_ncs_selection ,ncs_grad)
  # update newly created flex.vec3 with NCS copies
  for nrg in ncs_restraints_group_list:
    ncs_selection = nrg.master_iselection
    for ncs_copy in nrg.copies:
      asu_selection = ncs_copy.copy_iselection
      ncs_grad_portion = ncs_grad.select(ncs_selection)
      g.set_selected(asu_selection,ncs_grad_portion)
  return g.as_double()

def restraints_target_and_grads(
      restraints_manager,
      xray_structure,
      ncs_restraints_group_list,
      refine_selection,
      lbfgs_self,
      grad,
      refine_sites=False,
      refine_u_iso=False,
      refine_transformations=False,
      iso_restraints=None,
      use_hd=None):
  assert [refine_sites, refine_u_iso, refine_transformations].count(True)==1
  ef_grad = None
  ef = None
  if(restraints_manager is not None):
    if(refine_sites):
      ef = restraints_manager.energies_sites(
        sites_cart        = xray_structure.sites_cart(),
        compute_gradients = True)
      ef_grad = ef.gradients
    elif(refine_u_iso):
      ef = restraints_manager.energies_adp_iso(
        xray_structure    = xray_structure,
        parameters        = iso_restraints,
        use_u_local_only  = iso_restraints.use_u_local_only,
        use_hd            = use_hd,
        compute_gradients = True)
      ef_grad = ef.gradients
    elif(refine_transformations):
      ef = restraints_manager.energies_sites(
        sites_cart        = xray_structure.sites_cart(),
        compute_gradients = True)
      ef_grad = nu.compute_transform_grad(
        grad_wrt_xyz      = ef.gradients.as_double(),
        ncs_restraints_group_list = ncs_restraints_group_list,
        xyz_asu           = lbfgs_self.xray_structure.sites_cart(),
        x                 = lbfgs_self.x)
  if(ef is not None): return ef.target, ef_grad
  elif not grad: return 0, 0
  elif isinstance(grad,flex.vec3_double):
    return 0, flex.vec3_double([(0,0,0),]*len(grad))
  elif isinstance(grad,flex.double):
    return 0, flex.double((0,0,0)*(len(grad)//3))
  else: return None,None

class target_function_and_grads_real_space(object):
  """
  Real-space target and gradients evaluator
  """
  def __init__(
        self,
        map_data,
        xray_structure,
        ncs_restraints_group_list,
        refine_selection,
        real_space_gradients_delta,
        restraints_manager=None,
        data_weight=None,
        refine_sites=False,
        refine_transformations=False):
    adopt_init_args(self, locals())
    self.extended_ncs_selection = nu.get_extended_ncs_selection(
      ncs_restraints_group_list=ncs_restraints_group_list,
      refine_selection=refine_selection)

    self.unit_cell = self.xray_structure.unit_cell()
    # get selection to refine
    asu_size = xray_structure.scatterers().size()
    if refine_sites:
      # Use all atoms to refine
      self.selection = flex.bool(asu_size, refine_selection)
    elif refine_transformations:
      # use only NCS related atoms (Without the Master NCS copy)
      self.selection = nu.get_ncs_related_selection(
        ncs_restraints_group_list=ncs_restraints_group_list,
        asu_size=asu_size)

  def data_target_and_grads(self, compute_gradients, lbfgs_self):
    g = None
    tg = cctbx.maptbx.real_space_refinement_simple.target_and_gradients(
      unit_cell                  = self.unit_cell,
      density_map                = self.map_data,
      sites_cart                 = self.xray_structure.sites_cart(),
      real_space_gradients_delta = self.real_space_gradients_delta,
      selection                  = self.selection)
    t = tg.target()
    if compute_gradients:
      if self.refine_sites:
        g = tg.gradients()
      elif self.refine_transformations:
        grad_wrt_xyz = tg.gradients().as_double()
        g = nu.compute_transform_grad(
          grad_wrt_xyz      = grad_wrt_xyz,
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          xyz_asu           = lbfgs_self.xray_structure.sites_cart(),
          x                 = lbfgs_self.x)
    return t, g

  def target_and_gradients(self, compute_gradients, lbfgs_self, xray_structure):
    self.xray_structure.set_sites_cart(sites_cart = xray_structure.sites_cart())
    g = None
    t_data, g_data = self.data_target_and_grads(compute_gradients, lbfgs_self)
    t_restraints, g_restraints = restraints_target_and_grads(
      restraints_manager     = self.restraints_manager,
      xray_structure         = self.xray_structure,
      ncs_restraints_group_list = self.ncs_restraints_group_list,
      refine_selection       = self.refine_selection,
      refine_sites           = self.refine_sites,
      refine_transformations = self.refine_transformations,
      lbfgs_self             = lbfgs_self,
      grad                   = g_data)
    if(self.data_weight is None): self.data_weight=1.
    t = t_data*self.data_weight + t_restraints
    if(compute_gradients):
      g = g_data*self.data_weight + g_restraints
      if(not self.refine_transformations):
        g = grads_asu_to_one_ncs(
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          extended_ncs_selection    = self.extended_ncs_selection,
          grad                      = g,
          refine_sites              = self.refine_sites).as_double()
    return t, g

class target_function_and_grads_reciprocal_space(object):
  """
  Reciprocal-space target and gradients evaluator
  """
  def __init__(
        self,
        fmodel,
        ncs_restraints_group_list,
        refine_selection,
        restraints_manager=None,
        data_weight=None,
        refine_sites=False,
        refine_u_iso=False,
        refine_transformations=False,
        iso_restraints = None,
        use_hd         = False):
    adopt_init_args(self, locals())
    self.extended_ncs_selection = nu.get_extended_ncs_selection(
      ncs_restraints_group_list=ncs_restraints_group_list,
      refine_selection=refine_selection)
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.xray_structure = self.fmodel.xray_structure
    if self.refine_sites:
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)
    elif self.refine_u_iso:
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        u_iso      = True)
    elif self.refine_transformations:
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)

  def data_target_and_grads(self, compute_gradients, lbfgs_self):
    g = None
    tgx = self.x_target_functor(compute_gradients=compute_gradients)
    t = tgx.target_work()
    if compute_gradients:
      if self.refine_sites:
        g = flex.vec3_double(
          tgx.gradients_wrt_atomic_parameters(site=True).packed())
      elif self.refine_u_iso:
        g = tgx.gradients_wrt_atomic_parameters(u_iso=True)
      elif self.refine_transformations:
        grad_wrt_xyz = tgx.gradients_wrt_atomic_parameters(site=True).packed()
        g = nu.compute_transform_grad(
          grad_wrt_xyz      = grad_wrt_xyz,
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          xyz_asu           = lbfgs_self.xray_structure.sites_cart(),
          x                 = lbfgs_self.x)
    return t, g

  def target_and_gradients(self, compute_gradients, lbfgs_self, xray_structure):
    self.xray_structure.set_sites_cart(sites_cart = xray_structure.sites_cart())
    self.fmodel.update_xray_structure(
      xray_structure = self.xray_structure,
      update_f_calc  = True)
    g = None
    t_data, g_data = self.data_target_and_grads(compute_gradients, lbfgs_self)
    t_restraints, g_restraints =  restraints_target_and_grads(
      restraints_manager     = self.restraints_manager,
      xray_structure         = self.xray_structure,
      ncs_restraints_group_list = self.ncs_restraints_group_list,
      refine_selection       = self.refine_selection,
      refine_sites           = self.refine_sites,
      refine_u_iso           = self.refine_u_iso,
      refine_transformations = self.refine_transformations,
      lbfgs_self             = lbfgs_self,
      iso_restraints         = self.iso_restraints,
      use_hd                 = self.use_hd,
      grad                   = g_data)
    if(self.data_weight is None): self.data_weight=1.
    t = t_data*self.data_weight + t_restraints
    if(compute_gradients):
      g = g_data*self.data_weight + g_restraints
      if(not self.refine_transformations):
        g = grads_asu_to_one_ncs(
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          extended_ncs_selection    = self.extended_ncs_selection,
          grad                      = g,
          refine_sites              = self.refine_sites).as_double()
    return t, g

  def finalize(self):
    self.fmodel.xray_structure.tidy_us()
    self.fmodel.xray_structure.apply_symmetry_sites()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True,
      update_f_mask  = True)

class lbfgs(object):
  def __init__(self,
        ncs_restraints_group_list,
        target_and_grads_object,
        xray_structure,
        refine_selection             = None,
        finite_grad_differences_test = False,
        finite_grad_difference_val   = 0,
        max_iterations               = 35,
        refine_sites                 = False,
        refine_u_iso                 = False,
        refine_transformations       = False,
        use_strict_ncs               = True):
    """
    NCS constrained ADP and coordinates refinement. Also refines NCS operators.
    """
    adopt_init_args(self, locals())
    self.extended_ncs_selection = nu.get_extended_ncs_selection(
      ncs_restraints_group_list=ncs_restraints_group_list,
      refine_selection=refine_selection)
    assert [self.refine_sites,
            self.refine_u_iso, self.refine_transformations].count(True) == 1
    self.total_asu_length = len(xray_structure.sites_cart())
    assert self.refine_selection.size() > 0
    traditional_convergence_test_eps = 1.0e-6
    if self.use_strict_ncs:
      xray_structure_one_ncs_copy = xray_structure.select(
        self.extended_ncs_selection)
    else:
      xray_structure_one_ncs_copy = xray_structure
    if self.refine_sites:
      self.x = xray_structure_one_ncs_copy.sites_cart().as_double()
    elif self.refine_u_iso:
      assert xray_structure_one_ncs_copy.scatterers().size() == \
        xray_structure_one_ncs_copy.use_u_iso().count(True)
      self.x = xray_structure_one_ncs_copy.extract_u_iso_or_u_equiv()
    elif self.refine_transformations:
      self.x = nu.concatenate_rot_tran(
        ncs_restraints_group_list=self.ncs_restraints_group_list)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations,
        traditional_convergence_test_eps=traditional_convergence_test_eps),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    if(getattr(self.target_and_grads_object, "finalize", None)):
      self.target_and_grads_object.finalize()

  def compute_functional_and_gradients(self, compute_gradients=True):
    t,g = self.target_and_grads_object.target_and_gradients(
      compute_gradients=compute_gradients, lbfgs_self=self,
      xray_structure = self.update_xray_structure())
    if(self.finite_grad_differences_test and compute_gradients):
      self.finite_difference_test(g)
    return t, g

  def update_xray_structure(self, x=None):
    """
    Update xray_structure with refined parameters, then update
    fmodel object with updated xray_structure.
    """
    if not x: x = self.x
    if self.refine_transformations:
      # update the ncs_restraint_groups transforms
      self.ncs_restraints_group_list = nu.update_rot_tran(
        x=x, ncs_restraints_group_list=self.ncs_restraints_group_list)
      # Use the new transformations to create the ASU
      x_ncs = nu.get_ncs_sites_cart(self).as_double()
      x_asu = self.refinable_params_one_ncs_to_asu(x_ncs)
      self.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(x_asu))
    elif self.refine_sites:
      x_asu = self.refinable_params_one_ncs_to_asu(x)
      self.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(x_asu))
    elif self.refine_u_iso:
      x_asu = self.refinable_params_one_ncs_to_asu()
      self.xray_structure.set_u_iso(values = x_asu)
    return self.xray_structure

  def refinable_params_one_ncs_to_asu(self, x=None):
    """
    Expand refinabale parameters corresponding to one NCS copy to parameters
    corresponding to whole ASU.
    """
    x_old = self.xray_structure.sites_cart()
    if not x : x = self.x
    if self.refine_sites or self.refine_transformations:
      if self.use_strict_ncs or self.refine_transformations:
        new_x = nu.apply_transforms(
          ncs_coordinates = flex.vec3_double(x),
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          total_asu_length = x_old.size(),
          extended_ncs_selection = self.extended_ncs_selection,
          round_coordinates = False)
        new_x = new_x.select(self.refine_selection)
        new_x = x_old.set_selected(self.refine_selection,new_x)
      return new_x.as_double()
    elif self.refine_u_iso:
      if self.use_strict_ncs:
        return grads_one_ncs_to_asu(
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          extended_ncs_selection = self.extended_ncs_selection,
          total_asu_length = x_old.size(),
          ncs_grad = x)
      else:
        return flex.double(list(x))

  def finite_difference_test(self,g):
    """
    Compare analytical and finite differences gradients.

    finite_grad_difference_val = abs(analytical - finite differences)
    """
    g = g.as_double()
    # find the index of the max gradient value
    # i_g_max = flex.max_index(flex.abs(g))

    # Set displacement for finite gradient calculation to 1e-5 of largest value
    if self.refine_transformations:
      d = 1.0e-3
      i_g_max = 0
    else:
      i_g_max = flex.max_index(flex.abs(g))
      d = max(self.x[i_g_max]*1e-6,1e-6)
    # calc t(x+d)
    self.x[i_g_max] = self.x[i_g_max] + d
    t1,_ = self.compute_functional_and_gradients(compute_gradients=False)
    # calc t(x-d)
    self.x[i_g_max] = self.x[i_g_max] - 2*d
    t2,_ = self.compute_functional_and_gradients(compute_gradients=False)
    # Return fmodel to the correct coordinates values
    self.x[i_g_max] = self.x[i_g_max] + d
    self.update_xray_structure()
    finite_gard = (t1-t2)/(d*2)
    self.finite_grad_difference_val = abs(g[i_g_max] - finite_gard)

  def group_same_transforms(self):
    """
    Analyze and combine similar transformations
    """
    # todo: check if still needed
    from  iotbx.pdb.multimer_reconstruction import ncs_group_object
    new_transforms_obj =  ncs_group_object()
    pass

  def __call__(self):
    # TODO: Check if using this instead of "self" as a parameter works
    return self
