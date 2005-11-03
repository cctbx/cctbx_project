import cctbx.xray.structure_factors
import cctbx.xray.minimization
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import adptbx
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
from stdlib import math
import math, sys
from mmtbx.refinement import print_statistics
from libtbx.test_utils import approx_equal


class lbfgs(object):

  def __init__(self, xray_gradient_flags,
                     xray_structure,
                     restraints_manager,
                     wx,
                     wc,
                     wu,
                     fmodel,
                     wilson_b,
                     lbfgs_termination_params=None,
                     cos_sin_table=True,
                     use_fortran=False,
                     verbose=0,
                     iso_restraints=None,
                     alpha_w=None,
                     beta_w=None):
    adopt_init_args(self, locals())
    if(self.fmodel.alpha_beta_params.method == "calc"):
       if(self.fmodel.alpha_beta_params.fix_scale_for_calc_option == None):
          self.scale_ml = self.fmodel.scale_ml()
       else:
          self.scale_ml = self.fmodel.alpha_beta_params.fix_scale_for_calc_option
    if(self.fmodel.alpha_beta_params.method == "est"):
       self.scale_ml = 1.0
    assert approx_equal(       self.xray_structure.sites_cart(),
                        self.fmodel.xray_structure.sites_cart())
    assert approx_equal(       self.xray_structure.scatterers().extract_u_iso(),
                        self.fmodel.xray_structure.scatterers().extract_u_iso())
    self.f_obs_w = self.fmodel.f_obs_w()
    self.f_obs = self.fmodel.f_obs
    self.xray_structure.tidy_us(u_min = 1.e-6)
    assert self.fmodel.sf_algorithm in ("fft","direct")
    self.target_name = self.fmodel.target_name
    assert self.target_name in ("ml","mlhl") or self.target_name.count("ls") == 1
    self.fmodel_copy = self.fmodel.deep_copy()
    if(self.alpha_w is None or self.beta_w is None):
       self.alpha_w, self.beta_w = self.fmodel_copy.alpha_beta_w()
    else:
       assert self.alpha_w.data().size() == self.f_obs_w.data().size()
       assert self.beta_w.data().size() == self.f_obs_w.data().size()
    self.structure_factor_gradients = cctbx.xray.structure_factors.gradients(
                                                 miller_set    = self.f_obs_w,
                                                 cos_sin_table = cos_sin_table)
    self.x = flex.double(xray_structure.n_parameters(xray_gradient_flags), 0)
    self._scatterers_start = xray_structure.scatterers()
    self._scattering_dict = xray_structure.scattering_dict()
    self._d_min = self.f_obs_w.d_min()
    self.first_target_value = None
    self._lock_for_line_search = False
    self.minimizer = scitbx.lbfgs.run(
      target_evaluator          = self,
      termination_params        = lbfgs_termination_params,
      use_fortran               = use_fortran,
      exception_handling_params = scitbx.lbfgs.exception_handling_parameters(
                         ignore_line_search_failed_step_at_lower_bound = True))
    self.apply_shifts()
    del self._scatterers_start
    del self._scattering_dict
    del self._d_min
    self._lock_for_line_search = False
    self.compute_target(compute_gradients = False, mean_displacements = None)
    del self._lock_for_line_search
    self.final_target_value = self.f

  def apply_shifts(self):
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      gradient_flags = self.xray_gradient_flags,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      scatterers_shifted[i_seq].site = crystal.correct_special_position(
        crystal_symmetry=self.xray_structure,
        special_op=site_symmetry_table.get(i_seq).special_op(),
        site_frac=scatterers_shifted[i_seq].site)
    self.xray_structure.replace_scatterers(scatterers=scatterers_shifted)
    return apply_shifts_result.mean_displacements

  def compute_target(self, compute_gradients, mean_displacements):
    self.stereochemistry_residuals = None
    self.fmodel_copy.update_xray_structure(self.xray_structure,
                                           update_f_calc            = True,
                                           update_f_mask            = False,
                                           update_f_ordered_solvent = False)
    if(self.target_name in ("ml","mlhl")):
       xrtfr = self.fmodel_copy.xray_target_functor_result(
                                             compute_gradients = True,
                                             alpha             = self.alpha_w,
                                             beta              = self.beta_w,
                                             scale_ml          = self.scale_ml,
                                             flag              = "work")
    if(self.target_name.count("ls") == 1):
       xrtfr = self.fmodel_copy.xray_target_functor_result(
                                                    compute_gradients = True,
                                                    flag              = "work")
    self.f = xrtfr.target() * self.wx
    if(compute_gradients):
       sf = self.structure_factor_gradients(
                        xray_structure     = self.xray_structure,
                        mean_displacements = mean_displacements,
                        miller_set         = self.f_obs_w,
                        d_target_d_f_calc  = xrtfr.derivatives(),
                        gradient_flags     = self.xray_gradient_flags,
                        n_parameters       = self.x.size(),
                        algorithm          = self.fmodel.sf_algorithm).packed()
       self.g = sf * self.wx
    if(self.xray_gradient_flags.site
         and self.restraints_manager is not None
         and self.wc > 0.0):
       self.stereochemistry_residuals = self.restraints_manager.energies_sites(
                       sites_cart        = self.xray_structure.sites_cart(),
                       compute_gradients = compute_gradients,
                       lock_for_line_search = self._lock_for_line_search)
       self._lock_for_line_search = True
       self.f += self.stereochemistry_residuals.target * self.wc
       if(compute_gradients):
          xray.minimization.add_gradients(
                scatterers     = self.xray_structure.scatterers(),
                gradient_flags = self.xray_gradient_flags,
                xray_gradients = self.g,
                site_gradients = self.stereochemistry_residuals.gradients * self.wc)
    if(self.xray_gradient_flags.u_iso
          and self.restraints_manager.geometry is not None
          and self.wu > 0.0
          and self.iso_restraints is not None):
       energies_adp_iso = self.restraints_manager.energies_adp_iso(
         xray_structure=self.xray_structure,
         parameters=self.iso_restraints,
         wilson_b=self.wilson_b,
         compute_gradients=compute_gradients)
       self.f += energies_adp_iso.target * self.wu
       if(compute_gradients):
          if(self.xray_gradient_flags.sqrt_u_iso):
             u_iso_gradients = 2*mean_displacements*energies_adp_iso.gradients
          elif(self.xray_gradient_flags.tan_b_iso_max != 0):
             u_iso_max = adptbx.b_as_u(self.xray_gradient_flags.tan_b_iso_max)
             u_iso_gradients = \
               (u_iso_max/math.pi/(flex.pow2(mean_displacements)+1.0)) * \
               energies_adp_iso.gradients
          else:
             u_iso_gradients = energies_adp_iso.gradients
          xray.minimization.add_gradients(
                            scatterers      = self.xray_structure.scatterers(),
                            gradient_flags  = self.xray_gradient_flags,
                            xray_gradients  = self.g,
                            u_iso_gradients = u_iso_gradients * self.wu)

  def callback_after_step(self, minimizer):
    self._lock_for_line_search = False
    if (self.verbose > 0):
      print "refinement.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()

  def compute_functional_and_gradients(self):
    mean_displacements = self.apply_shifts()
    self.compute_target(
      compute_gradients=True,
      mean_displacements=mean_displacements)
    if (self.first_target_value is None):
      self.first_target_value = self.f
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.f, self.g
