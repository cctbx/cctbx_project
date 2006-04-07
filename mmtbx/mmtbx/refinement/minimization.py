from cctbx import xray
from cctbx import crystal
from cctbx.xray.structure import structure as cctbx_xray_structure
from cctbx.array_family import flex
import scitbx.lbfgs
from libtbx import adopt_init_args
from stdlib import math
import sys, time
from libtbx.test_utils import approx_equal


class lbfgs(object):

  def __init__(self, restraints_manager,
                     fmodel,
                     wx,
                     wc,
                     lbfgs_termination_params=None,
                     cos_sin_table=True,
                     use_fortran=False,
                     verbose=0,
                     iso_restraints=None,
                     alpha_w=None,
                     beta_w=None):
    adopt_init_args(self, locals())
    self.xray_structure = self.fmodel.xray_structure
    xray.set_scatterer_grad_flags(scatterers = self.xray_structure.scatterers(),
                                  site       = True)
    self.echem_start = None
    self.exray_start = None
    self.echem_final = None
    self.exray_final = None
    self.r_work_start= self.fmodel.r_work()
    self.r_free_start= self.fmodel.r_free()
    if(self.fmodel.alpha_beta_params.method == "calc"):
       if(self.fmodel.alpha_beta_params.fix_scale_for_calc_option == None):
          self.scale_ml = self.fmodel.scale_ml()
       else:
          self.scale_ml = self.fmodel.alpha_beta_params.fix_scale_for_calc_option
    if(self.fmodel.alpha_beta_params.method == "est"):
       self.scale_ml = 1.0
    self.f_obs_w = self.fmodel.f_obs_w()
    self.xray_structure.tidy_us(u_min = 1.e-6)
    self.target_name = self.fmodel.target_name
    assert self.target_name in ("ml","mlhl") or self.target_name.count("ls") == 1
    if(self.target_name in ("ml","mlhl", "lsm")):
       if(self.alpha_w is None or self.beta_w is None):
          self.alpha_w, self.beta_w = self.fmodel.alpha_beta_w()
       else:
          assert self.alpha_w.data().size() == self.f_obs_w.data().size()
          assert self.beta_w.data().size() == self.f_obs_w.data().size()
    self.x = flex.double(self.xray_structure.n_parameters(), 0)
    self._scatterers_start = self.xray_structure.scatterers()
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
    self._lock_for_line_search = False
    self.compute_target(compute_gradients = False)
    del self._lock_for_line_search
    self.final_target_value = self.f


  def apply_shifts(self):
    apply_shifts_result = xray.ext.minimization_apply_shifts(
      unit_cell      = self.xray_structure.unit_cell(),
      scatterers     = self._scatterers_start,
      shifts         = self.x)
    scatterers_shifted = apply_shifts_result.shifted_scatterers
    site_symmetry_table = self.xray_structure.site_symmetry_table()
    for i_seq in site_symmetry_table.special_position_indices():
      scatterers_shifted[i_seq].site = crystal.correct_special_position(
                crystal_symmetry = self.xray_structure,
                special_op       = site_symmetry_table.get(i_seq).special_op(),
                site_frac        = scatterers_shifted[i_seq].site)
    self.xray_structure.replace_scatterers(scatterers = scatterers_shifted)

  def compute_target(self, compute_gradients):
    self.stereochemistry_residuals = None
    self.fmodel.update_xray_structure(xray_structure = self.xray_structure,
                                      update_f_calc  = True)
    self.exray_final = self.fmodel.target_w(alpha    = self.alpha_w,
                                            beta     = self.beta_w,
                                            scale_ml = self.scale_ml)
    if(self.exray_start is None): self.exray_start = self.exray_final
    self.f = self.exray_final * self.wx
    if(compute_gradients):
       sf = self.fmodel.gradient_wrt_atomic_parameters(
                                                  alpha = self.alpha_w,
                                                  beta  = self.beta_w).packed()
       self.g = sf * self.wx
    if(self.restraints_manager is not None and self.wc > 0.0):
       self.stereochemistry_residuals = self.restraints_manager.energies_sites(
                       sites_cart        = self.xray_structure.sites_cart(),
                       compute_gradients = compute_gradients,
                       lock_for_line_search = self._lock_for_line_search)
       self._lock_for_line_search = True
       self.echem_final = self.stereochemistry_residuals.target
       if(self.echem_start is None): self.echem_start = self.echem_final
       self.f += self.echem_final * self.wc
       if(compute_gradients):
          xray.minimization.add_gradients(
             scatterers     = self.xray_structure.scatterers(),
             xray_gradients = self.g,
             site_gradients = self.stereochemistry_residuals.gradients*self.wc)

  def callback_after_step(self, minimizer):
    self._lock_for_line_search = False
    if (self.verbose > 0):
      print "refinement.minimization step: f,iter,nfun:",
      print self.f,minimizer.iter(),minimizer.nfun()

  def compute_functional_and_gradients(self):
    self.apply_shifts()
    self.compute_target(
      compute_gradients=True)
    if (self.first_target_value is None):
      self.first_target_value = self.f
    if (self.verbose > 1):
      print "xray.minimization line search: f,rms(g):",
      print self.f, math.sqrt(flex.mean_sq(self.g))
    return self.f, self.g

  def show(self, text=None, out=None):
    if (out is None): out = sys.stdout
    etotal_start = self.first_target_value
    etotal_final = self.final_target_value
    line_len = len("| "+text+"|")
    fill_len = 80 - line_len-1
    print >> out, "| "+text+"-"*(fill_len)+"|"
    print >> out, "| start r-factor (work) = %6.4f      final r-factor "\
     "(work) = %6.4f          |"%(self.r_work_start, self.fmodel.r_work())
    print >> out, "| start r-factor (free) = %6.4f      final r-factor "\
     "(free) = %6.4f          |"%(self.r_free_start, self.fmodel.r_free())
    print >> out, "|"+"-"*77+"|"
    format = "| %s = %s * %s + %s * %s"
    et = self.nn(etotal_start)
    wx = self.nn(self.wx)
    ex = self.nn(self.exray_start)
    wc = self.nn(self.wc)
    ec = self.nn(self.echem_start)
    h = "| T_start "+" "*abs(len("T_start")-len(et))+\
        "= wxc"+" "*abs(len("= wxc")-len(wx)-3)+\
        "* Exray "+" "*abs(len("* Exray ")-len(ex)-3)+\
        "+ wc "+" "*abs(len("+ wc ")-len(wc)-3)+\
        "* Echem "+" "*abs(len("* Echem ")-len(ec)-3)
    print >> out, h + " "*(78-len(h))+"|"
    st = format % (et,wx,ex,wc,ec)
    print >> out, st + " "*(78-len(st))+"|"

    et = self.nn(etotal_final)
    wx = self.nn(self.wx)
    ex = self.nn(self.exray_final)
    wc = self.nn(self.wc)
    ec = self.nn(self.echem_final)
    print >> out, "|"+"-"*77+"|"
    h = "| T_final "+" "*abs(len("T_final")-len(et))+\
        "= wxc"+" "*abs(len("= wxc")-len(wx)-3)+\
        "* Exray "+" "*abs(len("* Exray ")-len(ex)-3)+\
        "+ wc "+" "*abs(len("+ wc ")-len(wc)-3)+\
        "* Echem "+" "*abs(len("* Echem ")-len(ec)-3)
    print >> out, h + " "*(78-len(h))+"|"
    st = format % (et,wx,ex,wc,ec)
    print >> out, st + " "*(78-len(st))+"|"
    print >> out, "|"+"-"*77+"|"
    print >> out, "| number of iterations = %4d    |    number of function "\
                  "evaluations = %4d   |"%(self.minimizer.iter(),
                  self.minimizer.nfun())
    print >> out, "|"+"-"*77+"|"
    out.flush()

  def nn(self, x):
    x = str(x)
    try: return x[:len(x[:x.index(".")])+7]
    except: return "undef"
