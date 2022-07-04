from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from scitbx.array_family import flex
from libtbx import adopt_init_args
import mmtbx.ncs.ncs_utils as nu
from cctbx import xray
import scitbx.lbfgs
from six.moves import zip

def grads_asu_to_one_ncs(
      ncs_restraints_group_list,
      extended_ncs_selection,
      grad,
      refine_sites):
  """
  Apply NCS constraints to reduce the gradients corresponding to whole ASU, to
  gradients corresponding to one NCS copy.

  Args:
    ncs_restraints_group_list: list of ncs_restraint_group objects
    grad: gradient of the complete asu
    refine_sites: (bool) Flag indicating sites refinement
    extended_ncs_selection (flex.siz_t): selection of all ncs groups master ncs
      selection and non ncs related portions that are being refined (exclude
      NCS copies)
  """
  for nrg in ncs_restraints_group_list:
    ncs_selection = nrg.master_iselection
    for ncs_copy in nrg.copies:
      asu_selection = ncs_copy.iselection
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
                         master_grad):
  """
  Expand average gradient of a single NCS to all ASU
  (only for u_iso refinement)

  Args:
    ncs_restraints_group_list: list of ncs_restraint_group objects
    total_asu_length (int): length of the complete ASU
    extended_ncs_selection (flex.size_t): selection of all ncs groups master ncs
      selection and non ncs related portions that are being refined (exclude
      NCS copies)
    master_grad: gradient of a single ncs copy (the master copy)
  """
  g_length = total_asu_length
  if isinstance(master_grad,flex.vec3_double):
    g = flex.vec3_double([(0.0,0.0,0.0)]*g_length)
  elif isinstance(master_grad,flex.double):
    g = flex.double([0.0]*g_length)
  else:
    raise TypeError('Non supported grad type')
  # update newly created flex.vec3 with master NCS info
  g.set_selected(extended_ncs_selection ,master_grad)
  # update newly created flex.vec3 with NCS copies
  for nrg in ncs_restraints_group_list:
    master_selection = nrg.master_iselection
    master_grad_portion = g.select(master_selection)
    for ncs_copy in nrg.copies:
      copy_selection = ncs_copy.iselection
      g.set_selected(copy_selection,master_grad_portion)
  return g.as_double()

def restraints_target_and_grads(
      ncs_restraints_group_list,
      restraints_manager,
      xray_structure,
      grad,
      x,
      refine_sites=False,
      refine_u_iso=False,
      refine_transformations=False,
      iso_restraints=None,
      use_hd=None):
  """
  Args:
    ncs_restraints_group_list: NCS operator list
    restraints_manager: (object)
    xray_structure: (object)
    grad: gradient without restraints
    x: refined parameter
    refine_sites: (bool) indicate refinement type
    refine_u_iso: (bool) indicate refinement type
    refine_transformations: (bool) indicate refinement type
    iso_restraints: (object) iso restraints parameters
    use_hd: (bool) Use hydrogen

  Returns:
    ef_target: (float) restraints target function value
    ef_grad: (flex.double or flex.vec3_double) restraints gradient
  """
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
        xyz_asu           = xray_structure.sites_cart(),
        x                 = x)
  if(ef is not None): return ef.target, ef_grad
  elif not grad: return 0, 0
  elif isinstance(grad,flex.vec3_double):
    return 0, flex.vec3_double([(0,0,0),]*len(grad))
  elif isinstance(grad,flex.double):
    return 0, flex.double((0,0,0)*(len(grad)//3))
  else: return None,None


class target_function_and_grads_geometry_minimization(object):
  """
  Target and gradients evaluator for geometry minimization (no data)
  """
  def __init__(
        self,
        xray_structure,
        ncs_restraints_group_list,
        refine_selection=None,
        use_ncs_constraints=True,
        restraints_manager=None,
        refine_sites=False,
        refine_transformations=False):
    adopt_init_args(self, locals())
    self.refine_selection = nu.get_refine_selection(
      refine_selection=self.refine_selection,
      number_of_atoms=self.xray_structure.sites_cart().size())
    self.extended_ncs_selection = ncs_restraints_group_list.get_extended_ncs_selection(
        refine_selection=self.refine_selection)
    self.unit_cell = self.xray_structure.unit_cell()

  def target_and_gradients(self,compute_gradients,xray_structure,x):
    self.xray_structure.set_sites_cart(sites_cart = xray_structure.sites_cart())
    g = flex.vec3_double(xray_structure.scatterers().size(), [0,0,0])
    t = [0]*xray_structure.scatterers().size()
    t_restraints, g_restraints = restraints_target_and_grads(
      restraints_manager     = self.restraints_manager,
      xray_structure         = self.xray_structure,
      ncs_restraints_group_list = self.ncs_restraints_group_list,
      refine_sites           = self.refine_sites,
      refine_transformations = self.refine_transformations,
      x                      = x,
      grad                   = g)
    t = t_restraints
    if(compute_gradients):
      g = g_restraints
      if not self.refine_transformations:
        if self.use_ncs_constraints:
          g = grads_asu_to_one_ncs(
            ncs_restraints_group_list = self.ncs_restraints_group_list,
            extended_ncs_selection    = self.extended_ncs_selection,
            grad                      = g,
            refine_sites              = self.refine_sites).as_double()
        else:
          g = g.as_double()
    return t, g

class target_function_and_grads_real_space(object):
  """
  Real-space target and gradients evaluator
  """
  def __init__(
        self,
        map_data,
        xray_structure,
        ncs_restraints_group_list,
        real_space_gradients_delta,
        refine_selection=None,
        use_ncs_constraints=True,
        restraints_manager=None,
        data_weight=None,
        refine_sites=False):
    adopt_init_args(self, locals())
    self.refine_selection = nu.get_refine_selection(
      refine_selection=self.refine_selection,
      number_of_atoms=self.xray_structure.sites_cart().size())
    self.extended_ncs_selection = ncs_restraints_group_list.get_extended_ncs_selection(
        refine_selection=self.refine_selection)
    self.unit_cell = self.xray_structure.unit_cell()
    # get selection to refine
    asu_size = xray_structure.scatterers().size()
    self.selection = flex.bool(asu_size, refine_selection)

  def data_target_and_grads(self, compute_gradients,x):
    """
    Args:
      compute_gradients: (bool) when True compute gradients
      x: refined parameters
    """
    g = None
    tg = maptbx.target_and_gradients_simple(
      unit_cell   = self.unit_cell,
      map_target  = self.map_data,
      sites_cart  = self.xray_structure.sites_cart(),
      delta       = self.real_space_gradients_delta,
      selection   = self.selection)
    t = tg.target()*(-1)
    if compute_gradients:
      if self.refine_sites:
        g = tg.gradients()*(-1)
    return t, g

  def target_and_gradients(self,compute_gradients,xray_structure,x):
    """
    Args:
      compute_gradients: (bool) when True compute gradients
      xray_structure: xray_structure object
      x: refined parameters
    """
    self.xray_structure.set_sites_cart(sites_cart = xray_structure.sites_cart())
    g = None
    t_data, g_data = self.data_target_and_grads(compute_gradients,x)
    t_restraints, g_restraints = restraints_target_and_grads(
      restraints_manager     = self.restraints_manager,
      xray_structure         = self.xray_structure,
      ncs_restraints_group_list = self.ncs_restraints_group_list,
      refine_sites           = self.refine_sites,
      x                      = x,
      grad                   = g_data)
    if(self.data_weight is None): self.data_weight=1. #XXX BAD!
    t = t_data*self.data_weight + t_restraints
    if(compute_gradients):
      g = g_data*self.data_weight + g_restraints
      if self.use_ncs_constraints:
        g = grads_asu_to_one_ncs(
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          extended_ncs_selection    = self.extended_ncs_selection,
          grad                      = g,
          refine_sites              = self.refine_sites).as_double()
      else:
        g = g.as_double()
    return t, g

class target_function_and_grads_reciprocal_space(object):
  """
  Reciprocal-space target and gradients evaluator

  Args:
    refine_selection (flex.size_t): of all ncs related copies and
      non ncs related parts to be refined
  """
  def __init__(
        self,
        fmodel,
        ncs_restraints_group_list,
        refine_selection=None,
        use_ncs_constraints=True,
        restraints_manager=None,
        data_weight=None,
        refine_sites=False,
        refine_u_iso=False,
        refine_transformations=False,
        iso_restraints = None,
        use_hd         = False):
    adopt_init_args(self, locals())
    asu_size = self.fmodel.xray_structure.sites_cart().size()
    self.refine_selection = nu.get_refine_selection(
      refine_selection=self.refine_selection,
      number_of_atoms=asu_size)
    self.extended_ncs_selection = ncs_restraints_group_list.get_extended_ncs_selection(
        refine_selection=self.refine_selection)
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

  def data_target_and_grads(self, compute_gradients, sites_cart, x):
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
          xyz_asu           = sites_cart,
          x                 = x)
    return t, g

  def target_and_gradients(self, compute_gradients, xray_structure, x):
    sites_cart = xray_structure.sites_cart()
    self.xray_structure.set_sites_cart(sites_cart=sites_cart)
    self.fmodel.update_xray_structure(
      xray_structure = self.xray_structure,
      update_f_calc  = True)
    g = None
    t_data, g_data = self.data_target_and_grads(
      compute_gradients=compute_gradients, x=x, sites_cart=sites_cart)
    t_restraints, g_restraints =  restraints_target_and_grads(
      restraints_manager     = self.restraints_manager,
      xray_structure         = self.xray_structure,
      ncs_restraints_group_list = self.ncs_restraints_group_list,
      refine_sites           = self.refine_sites,
      refine_u_iso           = self.refine_u_iso,
      refine_transformations = self.refine_transformations,
      x                      = x,
      iso_restraints         = self.iso_restraints,
      use_hd                 = self.use_hd,
      grad                   = g_data)
    if(self.data_weight is None): self.data_weight=1.
    if(self.refine_transformations): # no restraints
      t = t_data
    else:
      t = t_data*self.data_weight + t_restraints
    if(compute_gradients):
      if(self.refine_transformations): # no restraints
        g = g_data
      else:
        g = g_data*self.data_weight + g_restraints
      if(not self.refine_transformations):
        if self.use_ncs_constraints:
          g = grads_asu_to_one_ncs(
            ncs_restraints_group_list = self.ncs_restraints_group_list,
            extended_ncs_selection    = self.extended_ncs_selection,
            grad                      = g,
            refine_sites              = self.refine_sites).as_double()
        else:
          g = g.as_double()
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
        refine_transformations       = False):
    """
    NCS constrained ADP and coordinates refinement. Also refines NCS operators.
    """
    adopt_init_args(self, args=locals(),exclude=['ncs_restraints_group_list'])
    self.x_previous = None
    self.refine_selection = nu.get_refine_selection(
      refine_selection=self.refine_selection,
      number_of_atoms=self.xray_structure.sites_cart().size())
    self.use_ncs_constraints = target_and_grads_object.use_ncs_constraints
    self.ncs_restraints_group_list = ncs_restraints_group_list.deep_copy()
    self.ncs_groups_coordinates_centers = []
    self.extended_ncs_selection = self.ncs_restraints_group_list.get_extended_ncs_selection(
        refine_selection=self.refine_selection)
    assert [self.refine_sites,
            self.refine_u_iso, self.refine_transformations].count(True) == 1
    self.total_asu_length = len(xray_structure.sites_cart())
    traditional_convergence_test_eps = 1.0e-6
    if self.use_ncs_constraints:
      xray_structure_one_ncs_copy = xray_structure.select(
        self.extended_ncs_selection)
    else:
      xray_structure_one_ncs_copy = xray_structure.select(self.refine_selection)
    if self.refine_sites:
      self.x = xray_structure_one_ncs_copy.sites_cart().as_double()
    elif self.refine_u_iso:
      assert xray_structure_one_ncs_copy.scatterers().size() == \
        xray_structure_one_ncs_copy.use_u_iso().count(True)
      self.x = xray_structure_one_ncs_copy.extract_u_iso_or_u_equiv()
    elif self.refine_transformations:
      # move refinable parameters to coordinate center
      self.ncs_groups_coordinates_centers = self.ncs_restraints_group_list.get_ncs_groups_centers(
          sites_cart=self.xray_structure.sites_cart())
      self.ncs_restraints_group_list = self.ncs_restraints_group_list.shift_translation_to_center(
          shifts = self.ncs_groups_coordinates_centers)
      self.x = self.ncs_restraints_group_list.concatenate_rot_tran()
    lbfgs_core_params = scitbx.lbfgs.core_parameters(
      stpmax = 25.0)
    minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      core_params=lbfgs_core_params,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations,
        traditional_convergence_test_eps=traditional_convergence_test_eps),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_step_at_upper_bound=True,
        ignore_line_search_failed_maxfev=True))
    # change transforms to the original coordinate system
    if self.refine_transformations:
      self.ncs_restraints_group_list = self.ncs_restraints_group_list.shift_translation_back_to_place(
          shifts = self.ncs_groups_coordinates_centers)
    if(getattr(self.target_and_grads_object, "finalize", None)):
      self.target_and_grads_object.finalize()
    # pass the refined ncs_restraints_group_list to original object
    for g1,g2 in zip(ncs_restraints_group_list,self.ncs_restraints_group_list):
      for tr1,tr2 in zip(g1.copies,g2.copies):
        tr1.r = tr2.r
        tr1.t = tr2.t

  def compute_functional_and_gradients(self, compute_gradients=True):
    x_current = self.x
    if(self.x_previous is None):
      self.x_previous = x_current.deep_copy()
    else:
      xray.ext.damp_shifts(previous=self.x_previous, current=x_current,
        max_value=10.)
      self.x_previous = x_current.deep_copy()

    t,g = self.target_and_grads_object.target_and_gradients(
      compute_gradients=compute_gradients,
      xray_structure = self.update_xray_structure(),
      x=self.x)
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
      self.ncs_restraints_group_list.update_rot_tran(x=x)
      # Use the new transformations to create the ASU
      x_ncs = self.xray_structure.sites_cart().\
          select(self.extended_ncs_selection).as_double()
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
      if self.use_ncs_constraints or self.refine_transformations:
        new_x = nu.apply_transforms(
          ncs_coordinates = flex.vec3_double(x),
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          total_asu_length = x_old.size(),
          extended_ncs_selection = self.extended_ncs_selection,
          round_coordinates = False,
          center_of_coordinates = self.ncs_groups_coordinates_centers)
        new_x = new_x.select(self.refine_selection)
        new_x = x_old.set_selected(self.refine_selection,new_x)
      else:
        new_x = self.x
      return new_x.as_double()
    elif self.refine_u_iso:
      if self.use_ncs_constraints:
        return grads_one_ncs_to_asu(
          ncs_restraints_group_list = self.ncs_restraints_group_list,
          extended_ncs_selection = self.extended_ncs_selection,
          total_asu_length = x_old.size(),
          master_grad = x)
      else:
        return flex.double(list(x))

  def finite_difference_test(self,g):
    """
    Compare analytical and finite differences gradients.

    finite_grad_difference_val = abs(analytical - finite differences)
    """
    g = g.as_double()
    # find the index of the max gradient value
    i_g_max = flex.max_index(flex.abs(g))
    # Set displacement for finite gradient calculation
    d = max(self.x[i_g_max]*1e-6,1e-6)
    # calc t(x+d)
    self.x[i_g_max] += d
    t1,_ = self.compute_functional_and_gradients(compute_gradients=False)
    # calc t(x-d)
    self.x[i_g_max] -= 2*d
    t2,_ = self.compute_functional_and_gradients(compute_gradients=False)
    # Return fmodel to the correct coordinates values
    self.x[i_g_max] += d
    self.update_xray_structure()
    finite_gard = (t1-t2)/(d*2)
    self.finite_grad_difference_val = abs(g[i_g_max] - finite_gard)
