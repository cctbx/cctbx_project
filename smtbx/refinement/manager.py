from cctbx import xray
from cctbx import sgtbx
from cctbx.array_family import flex
import scitbx.lbfgs
from smtbx.refinement.minimization import lbfgs
import math
from cctbx.eltbx import sasaki
from cctbx import adptbx

class manager(object):
  def __init__(self,
               f_obs = None,
               f_sq_obs = None,
               xray_structure = None,
               lambda_= None,
               max_cycles=50,
               verbose=1,
               log=None):
    assert [f_obs,f_sq_obs].count(None) == 1
    assert lambda_ is not None
    assert xray_structure is not None
    if f_obs:
      self.refinement_type = xray.amplitude
      f_obs.set_observation_type_xray_amplitude()
      self.f_obs = f_obs
      self.f_sq_obs = f_obs.f_as_f_sq()
    else:
      self.refinement_type = xray.intensity
      f_sq_obs.set_observation_type_xray_intensity()
      self.f_sq_obs = f_sq_obs
      self.f_obs = f_sq_obs.f_sq_as_f()
    assert self.f_sq_obs.is_real_array()
    self.xray_structure = xray_structure
    self.max_cycles = max_cycles
    self.verbose = verbose
    self.log = log

    for sc in self.xray_structure.scatterers():
      if sc.scattering_type in ('H','D'):continue
      fp_fdp = sasaki.table(sc.scattering_type).at_angstrom(lambda_)
      sc.fp = fp_fdp.fp()
      sc.fdp = fp_fdp.fdp()

  def start(self):
    """ Start the refinement """
    self.filter_reflections()
    #self.xs0 = self.xs
    self.set_refinement_flags()
    self.setup_refinement()
    self.start_refinement()

  def filter_reflections(self):
    f_sq_obs = self.f_sq_obs
    for i in xrange(f_sq_obs.size()):
      if f_sq_obs.data()[i] < -f_sq_obs.sigmas()[i]:
        f_sq_obs.data()[i] = -f_sq_obs.sigmas()[i]
    merging = f_sq_obs\
            .eliminate_sys_absent()\
            .merge_equivalents()
    f_sq_obs = merging.array()
    unique_reflections = f_sq_obs.size()
    if self.verbose:
      print "R(int) = %.4f   R(sigma) = %.4f" \
            %(merging.r_int(), merging.r_sigma())
      print "Total reflection: %i" %self.f_sq_obs.size()
      print "Unique reflections: %i" %unique_reflections
    f_obs = f_sq_obs.f_sq_as_f()
    self.f_sq_obs = f_sq_obs
    self.f_obs = f_obs

  def set_refinement_flags(self):
    for a in self.xray_structure.scatterers():
      a.flags.set_grad_site(True)
      if a.flags.use_u_iso() == True:
        a.flags.set_grad_u_iso(True)
        a.flags.set_grad_u_aniso(False)
      if a.flags.use_u_aniso()== True:
        a.flags.set_grad_u_aniso(True)
        a.flags.set_grad_u_iso(False)

  def setup_refinement(self):
    if self.refinement_type is xray.intensity:
      weighting = (xray.weighting_schemes.shelx_weighting(),None)[0]
      ls = xray.unified_least_squares_residual(self.f_sq_obs,
                                               weighting=weighting)
    else:
      ls = xray.unified_least_squares_residual(self.f_obs)
    self.ls = ls

  def start_refinement(self):
    minimisation = my_lbfgs(
      delegate=lambda mini: self.on_cycle_finished(self.xray_structure, mini),
      target_functor=self.ls,
      xray_structure=self.xray_structure,
      #structure_factor_algorithm="direct",
      cos_sin_table=True,
      lbfgs_termination_params = scitbx.lbfgs.termination_parameters(
        max_iterations=self.max_cycles),
      verbose=1,
      log=self.log,
    )
    self.minimisation = minimisation

  def on_cycle_finished(self, xs, minimiser):
    """ called after each iteration of the given minimiser, xs being
    the refined structure. """
    self.show_cycle_summary(minimiser)

  def show_cycle_summary(self, minimizer):
    if self.verbose:
      print "Refinement Cycle: %i" %(minimizer.iter())
      print "wR2: %.4f" %self.wR2()
      print "GooF = %.4f for %i parameters" %self.GooF(minimizer)

  def show_final_summary(self):
    """ only to be called after minimisation has finished. """
    if self.verbose:
      print "R1 = %.4f for %i Fo > 4sig(Fo)" %self.R1(),
      print "and %.4f for all %i data" %(self.R1(all_data=True))
      print "wR2 = %.4f" %self.wR2()
      print "GooF = %.4f for %i parameters" %self.GooF(self.minimisation.minimizer)

  def f_obs_minus_f_calc_map(self, resolution):
    f_obs=self.f_obs
    f_sq_obs = self.f_sq_obs
    f_sq_obs = f_sq_obs.eliminate_sys_absent().average_bijvoet_mates()
    f_obs = f_sq_obs.f_sq_as_f()
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_obs,
      cos_sin_table=True
    )
    f_calc = sf(self.xray_structure, f_obs).f_calc()
    fc2 = flex.norm(f_calc.data())
    fo2 = f_sq_obs.data()
    wfo2 = 1./flex.pow2(f_sq_obs.sigmas())
    K2 = flex.mean_weighted(fo2*fc2, wfo2)/flex.mean_weighted(fc2*fc2, wfo2)
    K2 = math.sqrt(K2)
    f_obs_minus_f_calc = f_obs.f_obs_minus_f_calc(1./K2, f_calc)
    return f_obs_minus_f_calc.fft_map(
      symmetry_flags=sgtbx.search_symmetry_flags(use_space_group_symmetry=False),
      resolution_factor=resolution,
    )

  def iter_scatterers(self):
    """ an iterator over tuples (label, xyz, u, u_eq, symbol) """
    for a in self.xray_structure.scatterers():
      label = a.label
      xyz = a.site
      symbol = a.scattering_type
      if a.flags.use_u_iso():
        u = (a.u_iso,)
        u_eq = u[0]
      if a.flags.use_u_aniso():
        u_cif = adptbx.u_star_as_u_cart(self.xray_structure.unit_cell(), a.u_star)
        u = u_cif
        u_eq = adptbx.u_star_as_u_iso(self.xray_structure.unit_cell(), a.u_star)

      yield label, xyz, u, u_eq, symbol

  def calculate_R1(self, f_obs):
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_obs,
      cos_sin_table=True
    )
    f_calc = sf(self.xray_structure, f_obs).f_calc()
    ls_function = xray.unified_least_squares_residual(f_obs)
    ls = ls_function(f_calc, compute_derivatives=False)
    k = ls.scale_factor()
    fc = flex.abs(f_calc.data())
    fo = f_obs.data()
    return flex.sum(flex.abs(k*fc - fo)) / flex.sum(fo)

  def R1(self, all_data=False):
    f_obs = self.f_obs
    #f_obs = f_obs.as_non_anomalous_array().merge_equivalents().array()
    if not all_data:
      strong = f_obs.data() > 4*f_obs.sigmas()
      f_obs = f_obs.select(strong)
    R1 = self.calculate_R1(f_obs)
    return R1, f_obs.size()

  def wR2(self):
    f_sq_obs = self.f_sq_obs
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_sq_obs,
      cos_sin_table=True
    )
    f_calc = sf(self.xray_structure, f_sq_obs).f_calc()
    ls_function = xray.unified_least_squares_residual(
      f_sq_obs,
      weighting=xray.weighting_schemes.shelx_weighting()
    )
    ls = ls_function(f_calc, compute_derivatives=False)
    weights = ls_function.weighting().weights
    k = ls.scale_factor()
    f_sq_calc = f_calc.norm()
    fc_sq = f_sq_calc.data()
    fo_sq = f_sq_obs.data()
    return math.sqrt(flex.sum(weights * flex.pow2(fo_sq - k * fc_sq)) /
                     flex.sum(weights * flex.pow2(fo_sq)))

  def GooF(self, minimizer):
    f_sq_obs = self.f_sq_obs
    sf = xray.structure_factors.from_scatterers(
      miller_set=f_sq_obs,
      cos_sin_table=True
    )
    f_calc = sf(self.xray_structure, f_sq_obs).f_calc()
    ls_function = xray.unified_least_squares_residual(
      f_sq_obs,
      weighting=xray.weighting_schemes.shelx_weighting()
    )
    ls = ls_function(f_calc, compute_derivatives=False)
    weights = ls_function.weighting().weights
    k = ls.scale_factor()
    f_sq_calc = f_calc.norm()
    fc_sq = f_sq_calc.data()
    fo_sq = f_sq_obs.data()
    return math.sqrt(flex.sum(weights * flex.pow2(fo_sq - k * fc_sq)) / (fo_sq.size() - minimizer.n())), minimizer.n()

class my_lbfgs(lbfgs):
  def __init__(self, delegate, **kwds):
    self.callback_after_step = delegate
    super(my_lbfgs, self).__init__(**kwds)
