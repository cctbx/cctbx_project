from __future__ import division
from cctbx import xray
import scitbx.lbfgs
from libtbx import adopt_init_args
from scitbx.array_family import flex

class lbfgs(object):
  def __init__(self,
        fmodel,
        rotation_matrices,
        translation_vectors,
        geometry_restraints_manager = None,
        data_weight = None,
        ncs_atom_selection = None,
        finite_grad_differences_test = False,
        finite_grad_difference_val = 0,
        max_iterations = 35,
        refine_sites = False,
        refine_u_iso = False,
        use_strict_ncs = True,
        iso_restraints = None,
        use_hd = False):
    """
    NCS constrained ADP and coordinates refinement.
    """
    adopt_init_args(self, locals())
    assert [self.geometry_restraints_manager,
            self.data_weight].count(None) in [0,2]
    assert [len(rotation_matrices) == len(translation_vectors)]
    assert [self.refine_sites, self.refine_u_iso].count(True) == 1
    assert ncs_atom_selection.count(True) > 0
    self.number_of_ncs_copies = len(rotation_matrices)+1
    self.fmodel.xray_structure.scatterers().flags_set_grads(state=False)
    self.x_target_functor = self.fmodel.target_functor()
    # xray structure of NCS chains for self.x
    if self.use_strict_ncs:
      xray_structure_one_ncs_copy = self.fmodel.xray_structure.select(
        ncs_atom_selection)
    else:
      xray_structure_one_ncs_copy = self.fmodel.xray_structure
    if(self.refine_sites):
      self.x = xray_structure_one_ncs_copy.sites_cart().as_double()
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        site       = True)
    if(self.refine_u_iso):
      assert xray_structure_one_ncs_copy.scatterers().size() == \
        xray_structure_one_ncs_copy.use_u_iso().count(True)
      self.x = xray_structure_one_ncs_copy.extract_u_iso_or_u_equiv()
      xray.set_scatterer_grad_flags(
        scatterers = self.fmodel.xray_structure.scatterers(),
        u_iso      = True)
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=scitbx.lbfgs.termination_parameters(
        max_iterations=max_iterations),
      exception_handling_params=scitbx.lbfgs.exception_handling_parameters(
        ignore_line_search_failed_rounding_errors=True,
        ignore_line_search_failed_step_at_lower_bound=True,
        ignore_line_search_failed_maxfev=True))
    self.fmodel.xray_structure.tidy_us()
    self.fmodel.xray_structure.apply_symmetry_sites()
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def compute_functional_and_gradients(self, compute_gradients=True):
    self.update_fmodel()
    g = None
    tgx = self.x_target_functor(compute_gradients=compute_gradients)
    t = tgx.target_work()
    if(self.refine_sites):
      if(self.geometry_restraints_manager is not None):
        es = self.geometry_restraints_manager.energies_sites(
          sites_cart        = self.fmodel.xray_structure.sites_cart(),
          compute_gradients = True)
        t = t*self.data_weight + es.target
      if(compute_gradients):
        gx = flex.vec3_double(
          tgx.gradients_wrt_atomic_parameters(site=True).packed())
        if(self.geometry_restraints_manager is not None):
          gx = gx*self.data_weight + es.gradients
        g = self.grads_asu_to_one_ncs(grad=gx).as_double()
    if(self.refine_u_iso):
      if(self.geometry_restraints_manager is not None):
        eadp = self.geometry_restraints_manager.energies_adp_iso(
          xray_structure = self.fmodel.xray_structure,
          parameters = self.iso_restraints,
          use_u_local_only = self.iso_restraints.use_u_local_only,
          use_hd = self.use_hd,
          compute_gradients = True)
        t = t*self.data_weight + eadp.target
      if(compute_gradients):
        gx = tgx.gradients_wrt_atomic_parameters(u_iso=True)
        if(self.geometry_restraints_manager is not None):
          gx = gx*self.data_weight + eadp.gradients
        g = self.grads_asu_to_one_ncs(grad=gx).as_double()
    if(self.finite_grad_differences_test and compute_gradients):
      self.finite_difference_test(g)
    return t, g

  def update_fmodel(self, x=None):
    """
    Update xray_structure with updated refinable parameters, then update
    fmodel object with updated xray_structure.
    """
    if not x: x = self.x
    if(self.refine_sites):
      x_asu = self.refinable_params_one_ncs_to_asu(x)
      self.fmodel.xray_structure.set_sites_cart(
        sites_cart = flex.vec3_double(x_asu))
    if(self.refine_u_iso):
      x_asu = self.refinable_params_one_ncs_to_asu()
      self.fmodel.xray_structure.set_u_iso(values = x_asu)
    self.fmodel.update_xray_structure(
      xray_structure = self.fmodel.xray_structure,
      update_f_calc  = True)

  def refinable_params_one_ncs_to_asu(self, x=None):
    """
    Expand refinabale parameters corresponding to one NCS copy to parameters
    corresponding to whole ASU.
    """
    if(self.refine_sites):
      assert x is not None
      new_x = list(x)
      x = flex.vec3_double(x)
      if self.use_strict_ncs:
        for r,t in zip(self.rotation_matrices, self.translation_vectors):
          tmp_x = r.elems*x + t
          new_x += list(tmp_x.as_double())
      return flex.double(new_x)
    if(self.refine_u_iso):
      if self.use_strict_ncs:
        return flex.double(list(self.x)*self.number_of_ncs_copies)
      else:
        return flex.double(list(self.x))

  def grads_asu_to_one_ncs(self, grad):
    """
    Apply NCS constraints to reduce gradients corresponding to whole ASU to
    gradients corresponding to one NCS copy.
    """
    if self.use_strict_ncs:
      # gradients of the first NCS copy
      ncs_end = len(grad)//self.number_of_ncs_copies
      assert ncs_end*self.number_of_ncs_copies==len(grad)
      g_ave = grad[:ncs_end]
      for i in range(self.number_of_ncs_copies-1):
        g = grad[ncs_end*(i+1):ncs_end*(i+2)]
        if(self.refine_sites):
          rt = self.rotation_matrices[i].transpose().elems
          g = rt*g
        g_ave += g
      g_ave = g_ave * (1./self.number_of_ncs_copies)
      if(self.refine_sites): g_ave = flex.vec3_double(g_ave)
      assert type(grad)==type(g_ave)
      return g_ave
    else:
      return grad

  def finite_difference_test(self,g):
    """
    Compare analytical and finite differences gradients.

    finite_grad_difference_val = abs(analytical - finite differences)
    """
    g = g.as_double()
    # find the index of the max gradient value
    i_g_max = flex.max_index(flex.abs(g))
    # Set displacement for finite gradient calculation to 1e-5 of largest value
    d = self.x[i_g_max]*0.00001
    # calc t(x+d)
    self.x[i_g_max] = self.x[i_g_max] + d
    t1,_ = self.compute_functional_and_gradients(compute_gradients=False)
    # calc t(x-d)
    self.x[i_g_max] = self.x[i_g_max] - 2*d
    t2,_ = self.compute_functional_and_gradients(compute_gradients=False)
    # Return fmodel to the correct coordinates values
    self.x[i_g_max] = self.x[i_g_max] + d
    self.update_fmodel()
    finite_gard = (t1-t2)/(d*2)
    self.finite_grad_difference_val = abs(g[i_g_max] - finite_gard)
