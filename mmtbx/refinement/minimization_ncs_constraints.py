from __future__ import division
import mmtbx.utils.ncs_utils as nu
from cctbx import xray
import scitbx.lbfgs
from libtbx import adopt_init_args
from scitbx.array_family import flex
import cctbx.maptbx.real_space_refinement_simple

def grads_asu_to_one_ncs(
      rotation_matrices,
      number_of_ncs_copies,
      grad,
      refine_sites):
  """
  Apply NCS constraints to reduce gradients corresponding to whole ASU, to
  gradients corresponding to one NCS copy.
  """
  # gradients of the first NCS copy
  ncs_end = len(grad)//number_of_ncs_copies
  assert ncs_end*number_of_ncs_copies==len(grad)
  g_ave = grad[:ncs_end]
  for i in range(number_of_ncs_copies-1):
    g = grad[ncs_end*(i+1):ncs_end*(i+2)]
    if(refine_sites):
      rt = rotation_matrices[i].transpose().elems
      g = rt*g
    g_ave += g
  g_ave = g_ave * (1./number_of_ncs_copies)
  if(refine_sites): g_ave = flex.vec3_double(g_ave)
  assert type(grad)==type(g_ave)
  return g_ave

def restraints_target_and_grads(
      restraints_manager,
      xray_structure,
      rotation_matrices,
      lbfgs_self,
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
        rotation_matrices = rotation_matrices,
        xyz_ncs           = nu.get_ncs_sites_cart(lbfgs_self),
        x                 = lbfgs_self.x)
  if(ef is not None): return ef.target, ef_grad
  else:               return None, None

class target_function_and_grads_real_space(object):
  """
  Real-space target and gradients evaluator
  """
  def __init__(
        self,
        map_data,
        xray_structure,
        rotation_matrices,
        translation_vectors,
        real_space_gradients_delta,
        restraints_manager=None,
        data_weight=None,
        refine_sites=False,
        refine_transformations=False):
    adopt_init_args(self, locals())
    assert [len(self.rotation_matrices) == len(self.translation_vectors)]
    self.number_of_ncs_copies = len(self.rotation_matrices)+1
    self.unit_cell = self.xray_structure.unit_cell()
    self.selection = flex.bool(xray_structure.scatterers().size(), True)

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
          rotation_matrices = self.rotation_matrices,
          xyz_ncs           = nu.get_ncs_sites_cart(lbfgs_self),
          x                 = lbfgs_self.x)
    return t, g

  def target_and_gradients(self, compute_gradients, lbfgs_self, xray_structure):
    self.xray_structure.set_sites_cart(sites_cart = xray_structure.sites_cart())
    g = None
    t_data, g_data = self.data_target_and_grads(compute_gradients, lbfgs_self)
    t_restraints, g_restraints = restraints_target_and_grads(
      restraints_manager     = self.restraints_manager,
      xray_structure         = self.xray_structure,
      rotation_matrices      = self.rotation_matrices,
      refine_sites           = self.refine_sites,
      refine_transformations = self.refine_transformations,
      lbfgs_self             = lbfgs_self)
    if(self.data_weight is None): self.data_weight=1.
    if([t_restraints, g_restraints].count(None)==0):
      t = t_data*self.data_weight + t_restraints
      if(compute_gradients):
        g = g_data*self.data_weight + g_restraints
        if(not self.refine_transformations):
          g = grads_asu_to_one_ncs(
            rotation_matrices    = self.rotation_matrices,
            number_of_ncs_copies = self.number_of_ncs_copies,
            grad                 = g,
            refine_sites         = self.refine_sites).as_double()
    else:
      t = t_data*self.data_weight
      if(compute_gradients):
        g = g_data*self.data_weight
        if(not self.refine_transformations):
          g = grads_asu_to_one_ncs(
            rotation_matrices    = self.rotation_matrices,
            number_of_ncs_copies = self.number_of_ncs_copies,
            grad                 = g,
            refine_sites         = self.refine_sites).as_double()
    return t, g

class target_function_and_grads_reciprocal_space(object):
  """
  Reciprocal-space target and gradients evaluator
  """
  def __init__(
        self,
        fmodel,
        rotation_matrices,
        translation_vectors,
        restraints_manager=None,
        data_weight=None,
        refine_sites=False,
        refine_u_iso=False,
        refine_transformations=False,
        iso_restraints = None,
        use_hd         = False):
    adopt_init_args(self, locals())
    assert [len(self.rotation_matrices) == len(self.translation_vectors)]
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    self.number_of_ncs_copies = len(self.rotation_matrices)+1
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
          rotation_matrices = self.rotation_matrices,
          xyz_ncs           = nu.get_ncs_sites_cart(lbfgs_self),
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
      rotation_matrices      = self.rotation_matrices,
      refine_sites           = self.refine_sites,
      refine_u_iso           = self.refine_u_iso,
      refine_transformations = self.refine_transformations,
      lbfgs_self             = lbfgs_self,
      iso_restraints         = self.iso_restraints,
      use_hd                 = self.use_hd)
    if(self.data_weight is None): self.data_weight=1.
    if([t_restraints, g_restraints].count(None)==0):
      t = t_data*self.data_weight + t_restraints
      if(compute_gradients):
        g = g_data*self.data_weight + g_restraints
        if(not self.refine_transformations):
          g = grads_asu_to_one_ncs(
            rotation_matrices    = self.rotation_matrices,
            number_of_ncs_copies = self.number_of_ncs_copies,
            grad                 = g,
            refine_sites         = self.refine_sites).as_double()
    else:
      t = t_data*self.data_weight
      if(compute_gradients):
        g = g_data*self.data_weight
        if(not self.refine_transformations):
          g = grads_asu_to_one_ncs(
            rotation_matrices    = self.rotation_matrices,
            number_of_ncs_copies = self.number_of_ncs_copies,
            grad                 = g,
            refine_sites         = self.refine_sites).as_double()
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
        rotation_matrices,
        translation_vectors,
        target_and_grads_object,
        xray_structure,
        ncs_atom_selection           = None,
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
    assert [len(rotation_matrices) == len(translation_vectors)]
    assert [self.refine_sites,
            self.refine_u_iso, self.refine_transformations].count(True) == 1
    assert ncs_atom_selection.count(True) > 0
    self.number_of_ncs_copies = len(rotation_matrices)+1
    traditional_convergence_test_eps = 1.0e-6
    if self.use_strict_ncs:
      xray_structure_one_ncs_copy = xray_structure.select(
        ncs_atom_selection)
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
        self.rotation_matrices,self.translation_vectors)
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
      self.rotation_matrices, self.translation_vectors\
        = nu.separate_rot_tran(x=x, ncs_copies=len(self.rotation_matrices))
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
    if self.refine_sites or self.refine_transformations:
      assert x is not None
      new_x = list(x)
      x = flex.vec3_double(x)
      if self.use_strict_ncs or self.refine_transformations:
        for r,t in zip(self.rotation_matrices, self.translation_vectors):
          tmp_x = r.elems*x + t
          new_x += list(tmp_x.as_double())
      return flex.double(new_x)
    elif self.refine_u_iso:
      if self.use_strict_ncs:
        return flex.double(list(self.x)*self.number_of_ncs_copies)
      else:
        return flex.double(list(self.x))

  def finite_difference_test(self,g):
    """
    Compare analytical and finite differences gradients.

    finite_grad_difference_val = abs(analytical - finite differences)
    """
    g = g.as_double()
    # find the index of the max gradient value
    i_g_max = flex.max_index(flex.abs(g))
    # Set displacement for finite gradient calculation to 1e-5 of largest value
    d = max(self.x[i_g_max]*1e-6,1e-8)
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
