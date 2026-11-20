from __future__ import absolute_import, division, print_function
import boost_adaptbx.boost.python as bp
ext = bp.import_ext("smtbx_refinement_least_squares_ext")
from smtbx_refinement_least_squares_ext import *


import smtbx.refinement.weighting_schemes # import dependency
from cctbx import xray
import libtbx.load_env
from libtbx import adopt_optional_init_args
from scitbx import linalg
from scitbx.lstbx import normal_eqns
from scitbx.array_family import flex
from smtbx.structure_factors import direct
from smtbx.refinement.restraints import origin_fixing_restraints
import math

# Above this much memory per worker, the packed accumulator is not offered to a
# threaded build however well it does on smaller ones. Each worker streams the
# whole of its private matrix once per reflection, so it is quick while that
# matrix stays in cache and turns over hard -- as a cliff, not a slope -- as
# soon as it does not. The buffered accumulator has no such boundary: it hands
# the work to a syrk, which the BLAS blocks for whatever cache it finds at
# runtime. So this is the line between an algorithm that adapts to the machine
# and one that cannot, and it has to be drawn conservatively.
#
# Expressed as bytes rather than as a parameter count, because bytes are what
# carries to another machine: the packed matrix is n(n+1)/2 doubles, and the
# question is whether that stays in cache.
#
# It is deliberately *not* read from the CPU that happens to be running. That
# would need a platform call apiece and would tune the refinement to one
# machine; and it does not need asking. What matters is last-level cache **per
# hardware thread**, since every worker holds one of these at once, and core
# counts and cache sizes have grown together, so that ratio has stayed within a
# narrow band across a long span of hardware -- roughly a megabyte a thread on
# an old quad-core, on a mid-range laptop and on a many-core desktop alike.
#
# The budget is set under the low end of that band, so it is safe on the oldest
# machine worth supporting rather than right on a modern one and wrong
# elsewhere. It gives up a little of what a large-cache machine could have
# used. The other end needs no judgement at all: a few thousand parameters is
# megabytes a worker and hundreds of megabytes across them, past the last-level
# cache of any machine that exists.

blas_2_parallel_max_matrix_bytes = 384 << 10


def blas_2_parallel_max_parameters(max_matrix_bytes=None):
  """ The largest n whose packed normal matrix fits the per-thread budget.

  n(n+1)/2 doubles <= budget, solved for n and rounded down. Arithmetic only:
  nothing here asks the operating system anything, and the answer is the same
  on every platform.
  """
  import math
  if max_matrix_bytes is None:
    max_matrix_bytes = blas_2_parallel_max_matrix_bytes
  doubles = max_matrix_bytes/8.0
  # n(n+1)/2 <= doubles  =>  n <= (sqrt(8 doubles + 1) - 1)/2
  return int((math.sqrt(8.0*doubles + 1.0) - 1.0)/2.0)


# Below this much work -- reflections x parameters -- sharing the reflection
# pass over threads costs more than it saves. There are two figures because
# there are two accumulators and they scale quite differently, and using one
# number for both is a mistake that shows up as a regression on exactly the
# structures which fall between them.
#
# The packed accumulator threads almost freely: each worker keeps a small matrix
# which stays in cache, so it pays off on anything but a trivial problem, and
# its crossover sits low.
#
# The buffered one does not. Its workers contend over the BLAS's own threading
# of the syrk, so it only starts paying at roughly ten times the work. Using one
# figure for both is a mistake that shows up as a regression on exactly the
# problems which fall between them.
blas_2_parallel_work_threshold = 2e5
blas_3_parallel_work_threshold = 2e6


def worth_parallelising(n_parameters, n_reflections, available_threads=None):
  """ Whether a build of this size should share its reflection pass out.

  Kept beside the accumulator choice deliberately: the two decisions are one
  decision. Which threshold applies depends on which accumulator will be used,
  and that depends on the parameter count, so asking the questions separately
  gets problems between the two thresholds wrong -- applying the low threshold
  to the buffered accumulator makes them markedly slower than serial.
  """
  if available_threads is not None and available_threads <= 1:
    return False
  work = float(n_reflections)*n_parameters
  if n_parameters <= blas_2_parallel_max_parameters():
    return work >= blas_2_parallel_work_threshold
  return work >= blas_3_parallel_work_threshold


def crystallographic_ls_class(non_linear_ls_with_separable_scale_factor=None,
                              n_parameters=None,
                              may_parallelise=None,
                              ml_target=None):
  """ Construct a class for crystallographic L.S. based on the given engine

  n_parameters and may_parallelise are the problem this class is about to be
  used on, where the caller knows them. They only ever pick between the two
  stock accumulators; an engine passed in explicitly is always honoured, and
  with neither given the choice is exactly what it has always been.

  ml_target selects the accumulator family. Maximum likelihood requires the
  fixed-scale one: the separable accumulator recovers its own scale by variable
  projection from yo.yc/yc.yc, whereas under a likelihood the effective
  observation already contains that scale, alpha and beta having been estimated
  with it folded in.
  """
  def blas_pair(fixed):
    """ The (level 3, level 2) accumulators of the required family """
    if fixed:
      return (getattr(normal_eqns,
                      "non_linear_ls_with_fixed_scale_factor_BLAS_3", None),
              normal_eqns.non_linear_ls_with_fixed_scale_factor_BLAS_2)
    return (getattr(normal_eqns,
                    "non_linear_ls_with_separable_scale_factor_BLAS_3", None),
            normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2)

  def get_base_class(non_linear_ls_with_separable_scale_factor):
    base_class = non_linear_ls_with_separable_scale_factor
    if not base_class:
      level_3, level_2 = blas_pair(fixed=bool(ml_target))
      try:
        from fast_linalg import env
        base_class = level_3 if (env.initialised and level_3) else level_2
      except Exception:
        base_class = level_2
      # A threaded build of a small enough problem is the one case where the
      # packed rank-1 update beats the buffered rank-k one. Each thread
      # accumulates into a matrix of its own, and while that matrix is small it
      # stays in cache, where BLAS 3's row buffer, its n x n scratch and the
      # per-thread syrks competing with each other all cost more than they save.
      if (may_parallelise
          and n_parameters is not None
          and n_parameters <= blas_2_parallel_max_parameters()):
        base_class = level_2
    #print("Chosen: " + str(base_class))
    return base_class

  class klass(get_base_class(non_linear_ls_with_separable_scale_factor)):

    non_linear_ls_engine = get_base_class(non_linear_ls_with_separable_scale_factor)

    default_weighting_scheme = mainstream_shelx_weighting
    weighting_scheme = "default"
    origin_fixing_restraints_type = (
      origin_fixing_restraints.atomic_number_weighting)
    f_mask = None
    restraints_manager=None
    n_restraints = None
    initial_scale_factor = None
    may_parallelise = False
    use_openmp = False
    max_memory = 300
    # Maximum-likelihood target: None for ordinary least squares, in which case
    # the settings below have no effect; "mlf" for the amplitude (Rice) target
    # and "mli" for the intensity one, convolved with the experimental error.
    ml_target = None
    # Test-set flags, a flex.bool true for the free reflections. Required:
    # alpha and beta estimated on the working set track the model as it
    # overfits, driving beta towards zero, whereupon the guard b <= 1e-3 in
    # mlf.h returns zero target and zero gradient for a whole shell silently.
    ml_free_flags = None
    ml_free_reflections_per_bin = 140
    # Estimate alpha and beta from the Fc the previous build already computed,
    # rather than computing it again first. See build_ml_data.
    ml_reuse_f_calc = True
    # Which scale the likelihood holds fixed within a cycle:
    #   'unit'     1, alpha having absorbed the Fo/Fc scale (the default)
    #   'given'    whatever the caller passed
    #   'optimal'  re-derived from the current Fc each cycle
    #
    # 'unit' is not a fallback, it is the convention mmtbx.max_lik and mlf.h
    # jointly expect. The estimator returns beta in *Fo* units, while mlf.h
    # forms b = beta*k^2; so any k but 1 rescales a variance that was never on
    # the Fc scale to begin with. Measured on jaca_A14: k = 1 refines
    # 0.0272 -> 0.0249 alongside least squares' 0.0244, while k = 0.22 - the
    # actual Fo/Fc scale, with alpha then coming out near 1 as theory says it
    # should - makes beta 20x too small, the weights 20x too large, and the
    # ADPs exceed the Debye-Waller limit within five cycles.
    #
    # So alpha absorbing the scale is correct here, and the tidier-looking
    # arrangement is the broken one. Kept selectable because that is a
    # statement about mmtbx's convention rather than about the mathematics,
    # and a different estimator would want a different answer.
    ml_scale_mode = 'unit'
    # Refmac's convention of folding the experimental variance into beta.
    # Off for the intensity target, which models that error explicitly.
    ml_add_sigma_squared_to_beta = None
    # How much the BLAS 3 accumulator may buffer before folding its rows into
    # the normal matrix. Zero means take it from max_memory, which is the
    # budget the whole build is meant to keep to.
    normal_matrix_buffer_bytes = 0

    @staticmethod
    def accumulator_buffer_bytes(n_parameters, max_memory_mb):
      """ What is left of the memory budget once the result is paid for.

      The accumulator holds the normal matrix twice over -- the full n x n the
      rank-k update writes into, and the packed copy it hands out -- and buffers
      rows in whatever remains. Chunking the rows costs a pass over the result
      per chunk, so the fewer chunks the better and the budget is worth
      spending; but it is a budget, and the result comes out of it first.

      Capped as well as floored. The buffer is reserved up front, so a generous
      budget would have it allocate the lot, which is the thing chunking is for
      avoiding. Past the cap there is nothing to buy anyway: the time is flat
      from a few hundred megabytes upwards.
      """
      if not max_memory_mb:
        return 0                      # no budget given: the accumulator decides
      n = int(n_parameters)
      result = (n*n + n*(n + 1)//2)*8
      return max(8 << 20,
                 min(512 << 20, int(max_memory_mb)*1048576 - result))

    def __init__(self, observations, reparametrisation,
                 one_h_linearisation=None, **kwds):
      # before adopt_optional_init_args, so these have to come out of kwds here
      buffer_bytes = kwds.get('normal_matrix_buffer_bytes',
                              klass.normal_matrix_buffer_bytes)
      if not buffer_bytes:
        buffer_bytes = klass.accumulator_buffer_bytes(
          reparametrisation.n_independents,
          kwds.get('max_memory', klass.max_memory))
      super(klass, self).__init__(
        reparametrisation.n_independents, True, buffer_bytes)
      self.observations = observations
      self.reparametrisation = reparametrisation
      adopt_optional_init_args(self, kwds)
      self.one_h_linearisation = one_h_linearisation
      if not self.one_h_linearisation:
        self.one_h_linearisation = f_calc_function_default(direct.f_calc_modulus_squared(
          self.xray_structure,
          disp_correction=reparametrisation.dispersion_radial))
      if self.weighting_scheme == "default":
        self.weighting_scheme = self.default_weighting_scheme()
      self.origin_fixing_restraint = self.origin_fixing_restraints_type(
        self.xray_structure.space_group())
      self.taken_step = None
      self.restraints_normalisation_factor = None

    @property
    def xray_structure(self):
      return self.reparametrisation.structure

    @property
    def twin_fractions(self):
      return self.reparametrisation.twin_fractions

    # A subclass which assembles the solvent contribution once and keeps it,
    # the mask being fixed for the length of a run, puts it here and mask_data()
    # hands that one out rather than rebuilding it.
    f_mask_data = None

    def mask_data(self):
      """ The solvent contribution, in the form the builders expect.

      Factored out of build_up so that anything else linearising the same
      problem computes the same structure factors, c.f. smtbx.refinement.cgls.
      Leaving it to each caller to assemble is how the two paths came to
      disagree: an omitted mask does not fail, it quietly refines a different
      structure.
      """
      if self.f_mask_data is not None:
        return self.f_mask_data
      if self.f_mask is None:
        return MaskData(flex.complex_double())
      return MaskData(self.observations, self.xray_structure.space_group(),
        self.observations.fo_sq.anomalous_flag(), self.f_mask.indices(),
        self.f_mask.data())

    def build_up(self, objective_only=False):
      f_mask_data = self.mask_data()

      fc_correction = self.reparametrisation.fc_correction
      if fc_correction is None:
        fc_correction = xray.dummy_fc_correction()

      def build_normal_eqns(scale_factor, weighting_scheme, objective_only,
                            ml_data=None):
        return ext.build_normal_equations(
          self,
          self.observations,
          f_mask_data,
          weighting_scheme,
          scale_factor,
          self.one_h_linearisation,
          self.reparametrisation.jacobian_transpose_matching_grad_fc(),
          fc_correction,
          objective_only,
          self.may_parallelise,
          self.use_openmp,
          self.max_memory,
          ml_data if ml_data is not None else ext.ml_data())

      bootstrap_f_calc = None
      if not self.finalised: #i.e. never been called
        self.reparametrisation.linearise()
        self.reparametrisation.store()
        scale_factor = self.initial_scale_factor
        if scale_factor is None: # we haven't got one from previous refinement
          result = build_normal_eqns(scale_factor=None,
                                     weighting_scheme=sigma_weighting(),
                                     objective_only=True)
          if self.ml_target is not None:
            # nothing has fitted one yet, and scale_factor() would answer with
            # its own default of one
            self.ml_crystallographic_scale = \
              self.crystallographic_scale_factor(result.observables())
          scale_factor = self.scale_factor()
          # this pass computed Fc for the current model, and alpha and beta
          # want exactly that: no reason for the likelihood to compute it again
          bootstrap_f_calc = self.observations.fo_sq.array(
            data=result.f_calc(), sigmas=None)
      else: # use scale factor from the previous step
        scale_factor = self.scale_factor()

      ml_data = None
      if self.ml_target is not None:
        self.check_maximum_likelihood_accumulator()
        ml_data = self.build_ml_data(f_mask_data, fc_correction, scale_factor,
                                     f_calc=bootstrap_f_calc)

      self.reset()
      result = build_normal_eqns(scale_factor,
                                 self.weighting_scheme,
                                 objective_only,
                                 ml_data)
      self.f_calc = self.observations.fo_sq.array(
        data=result.f_calc(), sigmas=None)
      self.fc_sq = self.observations.fo_sq.array(
        data=result.observables(), sigmas=None)\
          .set_observation_type_xray_intensity()
      if self.ml_target is not None:
        self.ml_crystallographic_scale = \
          self.crystallographic_scale_factor(self.fc_sq.data())
      self.weights = result.weights()
      self.objective_data_only = self.objective()
      self.chi_sq_data_only = self.chi_sq()
      if self.restraints_manager is not None:
        # Here we determine a normalisation factor to place the restraints on the
        # same scale as the average residual. This is the normalisation
        # factor suggested in Giacovazzo and similar to that used by shelxl.
        # (shelx manual, page 5-1).
        # The factor 2 comes from the fact that we minimize 1/2 sum w delta^2
        if self.restraints_normalisation_factor is None:
          self.restraints_normalisation_factor \
              = 2 * self.objective_data_only/(self.n_equations-self.n_parameters)
        linearised_eqns = self.restraints_manager.build_linearised_eqns(
          self.xray_structure, self.reparametrisation.parameter_map())
        jacobian = \
          self.reparametrisation.jacobian_transpose_matching(
            self.reparametrisation.mapping_to_grad_fc_all).transpose()
        self.reduced_problem().add_equations(
          linearised_eqns.deltas,
          linearised_eqns.design_matrix * jacobian,
          linearised_eqns.weights * self.restraints_normalisation_factor,
          optimise_for_tall_matrix=False)
        self.n_restraints = linearised_eqns.n_restraints()
        self.chi_sq_data_and_restraints = self.chi_sq()
      if not objective_only:
        self.origin_fixing_restraint.add_to(
          self.step_equations(),
          self.reparametrisation.jacobian_transpose_matching_grad_fc(),
          self.reparametrisation.asu_scatterer_parameters)
      # Whether what is now standing is a Gauss-Newton system or only the
      # objective. An objective-only build still accumulates restraints and
      # origin fixing, so "is anything in there" does not answer it and the
      # difference is invisible until a Cholesky fails on the first parameter
      # no restraint happens to touch. Recorded here because build_up is the
      # only place that knows.
      self.normal_equations_are_complete = not objective_only
      # An override says which scale factor to weight the data with when there
      # are no fresh normal equations to take one from. There now are, so it
      # has served its purpose, and anything reading scale_factor() after this
      # -- the journal in normal_eqns_solving, for one -- wants the new one.
      self.overridden_scale_factor = None

    def check_maximum_likelihood_accumulator(self):
      """ Refuse a separable-scale accumulator under a likelihood target

      The separable accumulator solves for the scale by variable projection.
      Under a likelihood the observation it projects is the effective one,
      which already contains that scale, so the projection would solve for a
      quantity that is not free and would do so without failing. Callers
      passing an engine explicitly bypass the selection in
      crystallographic_ls_class, making this the only point at which the
      mistake can be detected.
      """
      engine = getattr(self, "non_linear_ls_engine", None)
      name = getattr(engine, "__name__", "") or ""
      if "separable" in name:
        raise RuntimeError(
          "maximum likelihood (ml_target=%r) needs a fixed-scale accumulator, "
          "but this refinement was built on %s. Construct the class with "
          "crystallographic_ls_class(ml_target=...), which picks the right "
          "family and BLAS level, or pass "
          "non_linear_ls_with_fixed_scale_factor_BLAS_2/_3 explicitly."
          % (self.ml_target, name))

    def build_ml_data(self, f_mask_data, fc_correction, scale_factor,
                      f_calc=None):
      """ alpha, beta, epsilon and the centric flags for one macro-cycle

      Re-estimated each cycle rather than once, alpha and beta measuring the
      discrepancy between model and truth which refinement is itself altering.

      alpha and beta have to be known *before* the build that uses them, while
      Fc for the current model only exists *after* one, so something has to
      give. Estimating them from the Fc the last build already computed costs
      nothing and lags them by one cycle; computing Fc again first costs a
      complete extra pass over every reflection - measured at 19% of an ML
      refinement on crambin - and buys alpha and beta that describe the model
      exactly as it stands.

      The lag is what Refmac and Phenix live with, sigma_A being estimated from
      the model as it enters a cycle. It is also small where it matters: alpha
      and beta are smoothed over resolution shells of some hundred reflections,
      so one cycle of shifts moves them by far less than the binning already
      does, and less still as the refinement converges and the shifts shrink.
      Set ml_reuse_f_calc False to pay for the extra pass instead.
      """
      from smtbx.refinement import sigma_a
      if self.ml_free_flags is None:
        raise RuntimeError(
          "maximum likelihood needs a free set: estimating alpha and beta on "
          "the working set follows the model as it overfits, and beta then "
          "falls below the target's own floor, which silently zeroes both the "
          "target and its gradient for whole resolution shells")
      if f_calc is None and self.ml_reuse_f_calc:
        f_calc = getattr(self, 'f_calc', None)
        # a model whose scatterer count changed under us invalidates it
        if (f_calc is not None
            and f_calc.size() != self.observations.fo_sq.size()):
          f_calc = None
      if f_calc is None:
        # nothing to reuse: the first cycle of a run started from a stored
        # scale factor, or reuse switched off
        #
        # self is the accumulator and has been finalised - by the scale
        # bootstrap on the first cycle, by the previous cycle after that. The
        # probe accumulates into it too, so it has to start clean; the real
        # build resets again immediately afterwards.
        self.reset()
        probe = ext.build_normal_equations(
          self,
          self.observations,
          f_mask_data,
          sigma_weighting(),
          scale_factor,
          self.one_h_linearisation,
          self.reparametrisation.jacobian_transpose_matching_grad_fc(),
          fc_correction,
          True,               # objective only: no gradients wanted here
          self.may_parallelise,
          self.use_openmp,
          self.max_memory,
          ext.ml_data())
        f_calc = self.observations.fo_sq.array(
          data=probe.f_calc(), sigmas=None)
      f_obs = self.observations.fo_sq.f_sq_as_f()
      add_sigma_sq = self.ml_add_sigma_squared_to_beta
      if add_sigma_sq is None:
        # the intensity target convolves with the experimental error already
        add_sigma_sq = (self.ml_target != "mli")
      # The scale, and Fc put onto the scale of Fo before alpha is asked for.
      #
      # This ordering is the whole of it. mlf.h computes a = alpha*k, so alpha
      # is defined against a *scaled* Fc: <Fo> = alpha <k Fc>. Estimating it
      # against an unscaled one instead leaves alpha absorbing the scale, and
      # the target then applies that scale a second time.
      #
      # It showed as <alpha> = 0.25 on a deposited model, where anything near a
      # correct structure should give something close to 1 at low angle; 0.25
      # is about the reciprocal of this data's Fo/Fc scale, which is alpha
      # standing in for it. The refinement moved atoms by 2 A per cycle without
      # improving R1, because every effective observation fo/alpha was four
      # times what it should have been.
      if self.ml_scale_mode == 'unit':
        self.ml_scale_factor = 1.0
      elif self.ml_scale_mode == 'given':
        self.ml_scale_factor = scale_factor if scale_factor else 1.0
      else:
        # on amplitudes, alpha and beta being defined on those
        fo = f_obs.data()
        fc = flex.abs(f_calc.data())
        denom = flex.sum(fc*fc)
        self.ml_scale_factor = (flex.sum(fo*fc)/denom if denom > 0
                                else (scale_factor or 1.0))
      k = self.ml_scale_factor
      alpha, beta = sigma_a.alpha_beta(
        f_obs=f_obs,
        f_calc=f_calc.customized_copy(data=f_calc.data()*k),
        r_free_flags=self.ml_free_flags,
        free_reflections_per_bin=self.ml_free_reflections_per_bin,
        add_sigma_squared_to_beta=add_sigma_sq)
      centric, epsilons = sigma_a.centric_flags_and_epsilons(f_obs)
      self.ml_alpha, self.ml_beta = alpha, beta
      # The free set goes in so the build can hold it out of the target sum.
      # It is what alpha and beta were just estimated from, and refining
      # against it too would make that estimate circular.
      free = self.ml_free_flags
      return ext.ml_data(alpha=alpha, beta=beta, epsilon=epsilons,
                         centric=centric,
                         intensity=(self.ml_target == "mli"),
                         free_flags=(free if free is not None
                                     else flex.bool(alpha.size(), False)),
                         scale_factor=k)

    # The Fo^2/Fc^2 scale of the last build under a likelihood target; None
    # until one has happened. See scale_factor().
    ml_crystallographic_scale = None

    def crystallographic_scale_factor(self, fc_sq):
      """ The Fo^2/Fc^2 scale, fitted here rather than read off the accumulator.

      Maximum likelihood holds the accumulator's scale at one - alpha carries
      the scale, so projecting out a second one would not be solving for
      anything free - which leaves optimal_scale_factor() an algebraic constant
      rather than the scale of the structure. R1, wR2, the difference map, the
      standard uncertainties and the CIF all read scale_factor(), and one in
      place of the true 0.08 makes them report a destroyed structure that is
      not destroyed: crambin came back at R1 2.25 from a model whose real R1
      was 0.33.

      Same variable projection the separable accumulator does, with the same
      weighting scheme, so the number means what it means everywhere else. The
      SHELX weights depend on the scale through their P term, so the projection
      is repeated once from its own answer; a third pass moves it by nothing
      that shows in a reported figure.
      """
      fo_sq = self.observations.fo_sq
      yo = fo_sq.data()
      k = self.ml_crystallographic_scale or 1.0
      for _ in range(2):
        w = self.weighting_scheme(yo, fo_sq.sigmas(), fc_sq,
                                  fo_sq.indices(), k)
        denom = flex.sum(w*fc_sq*fc_sq)
        if denom <= 0:
          return 1.0
        k = flex.sum(w*yo*fc_sq)/denom
      return k

    def parameter_vector_norm(self):
      return self.reparametrisation.norm_of_independent_parameter_vector

    # Set by a solver which determines the scale factor without building the
    # normal equations, c.f. smtbx.refinement.cgls; None means take it from
    # them as usual. build_up clears it, having made it obsolete.
    overridden_scale_factor = None

    def scale_factor(self):
      if self.overridden_scale_factor is not None:
        return self.overridden_scale_factor
      if self.ml_target is not None:
        # the accumulator's own is fixed at one by construction
        return (self.ml_crystallographic_scale
                if self.ml_crystallographic_scale is not None else 1.0)
      return self.optimal_scale_factor()

    def apply_shifts(self, shifts):
      """ Move the structure by the given increment of the independent
          parameters.

      Factored out of step_forward so that minimisers which do not follow the
      step obtained from the normal equations can reuse it, c.f.
      scitbx.lstbx.scipy_iterations.
      """
      self.reparametrisation.apply_shifts(shifts)
      self.reparametrisation.linearise()
      self.reparametrisation.store()
      self.taken_step = shifts.deep_copy()

    def step_forward(self):
      self.apply_shifts(self.step())

    def parameter_vector(self):
      """ The independent parameters, as apply_shifts indexes them.

      apply_shifts does not always move the parameters by exactly the shifts
      it is given, validate() constraining some of them -- a U_iso or an
      occupancy driven negative, an extinction parameter, the thickness. A
      minimiser which places the parameters itself has to be able to see where
      they really went, or every subsequent shift is out by the difference.
      """
      return self.reparametrisation.independent_parameter_vector()

    def step_backward(self):
      self.reparametrisation.apply_shifts(-self.taken_step)
      self.reparametrisation.linearise()
      self.reparametrisation.store()
      self.taken_step = None

    def goof(self):
      if self.ml_target is not None:
        return self.least_squares_statistics()[1]
      return math.sqrt(self.chi_sq_data_only)

    def least_squares_statistics(self):
      """ wR2 and the goodness of fit as a least-squares run would report them.

      Under a likelihood the accumulator's objective is over effective
      observations weighted by a curvature, so a wR2 or a goodness of fit built
      from it is not comparable with a least-squares run: crambin reports 0.08
      where a fitted model should be near one, and the number moves when the
      target changes rather than when the model does.

      These are the ordinary quantities - measured intensities, the
      refinement's own weighting scheme, the crystallographic scale - so
      choosing between targets compares like with like. The likelihood is still
      what is being minimised; this only describes the model it reaches.
      """
      fo_sq = self.observations.fo_sq
      yo, yc = fo_sq.data(), self.fc_sq.data()
      k = self.scale_factor()
      w = self.weighting_scheme(yo, fo_sq.sigmas(), yc, fo_sq.indices(), k)
      r = yo - k*yc
      weighted = flex.sum(w*r*r)
      denom = flex.sum(w*yo*yo)
      wr2 = math.sqrt(weighted/denom) if denom > 0 else 0.
      dof = yo.size() - self.n_parameters
      return wr2, (math.sqrt(weighted/dof) if dof > 0 else 0.)

    def variance_goof_factor(self):
      """ What a variance from the inverse normal matrix has to be scaled by.

      Least squares knows its weights only up to a scale - the SHELX scheme is
      not 1/sigma^2 - so the inverse matrix is corrected by the goodness of fit,
      which measures exactly that discrepancy. This is the usual practice and
      what every reported s.u. rests on.

      A likelihood does not need it and is damaged by it. Its weight is the
      Gauss-Newton curvature of -logL, 2a^2/(epsilon b), in absolute units, so
      the inverse matrix is already an inverse information and the goodness of
      fit is not a correction to anything - it is a quantity with no reason to
      approach one, measured at 0.08 on crambin and 11.8 on the same structure
      without restraints. Applying it scaled every s.u. by between 1/140 and
      140 depending on which, which in turn is what ShelXL's DAMP divides the
      shifts by, so the step limiter was working in units that moved under it.

      The curvature retained is the positive part only, the term
      (2 a fo/epsilon b)^2 X'(u) having been dropped to keep the matrix
      positive definite, so these s.u. are still systematically small. That is a
      known bias of a fixed sign and not a moving scale.
      """
      if self.ml_target is not None:
        return 1.0
      return self.restrained_goof()**2

    def restrained_goof(self):
      if self.restraints_manager is None:
        return self.goof()
      return math.sqrt(self.chi_sq_data_and_restraints)

    def wR2(self, cutoff_factor=None):
      if cutoff_factor is None:
        if self.ml_target is not None:
          return self.least_squares_statistics()[0]
        return math.sqrt(2*self.objective_data_only)
      fo_sq = self.observations.fo_sq
      strong = fo_sq.data() >= cutoff_factor*fo_sq.sigmas()
      fo_sq = fo_sq.select(strong)
      fc_sq = self.fc_sq.select(strong)
      wght = self.weights.select(strong)
      fc_sq = fc_sq.data()
      fo_sq = fo_sq.data()
      fc_sq *= self.scale_factor()
      wR2 = flex.sum(wght*flex.pow2((fo_sq-fc_sq)))/flex.sum(wght*flex.pow2(fo_sq))
      return math.sqrt(wR2)

    def r1_factor(self, cutoff_factor=None):
      fo_sq = self.observations.fo_sq
      if cutoff_factor is not None:
        strong = fo_sq.data() >= cutoff_factor*fo_sq.sigmas()
        fo_sq = fo_sq.select(strong)
        fc_sq = self.fc_sq.select(strong)
      else:
        fc_sq = self.fc_sq
      f_obs = fo_sq.f_sq_as_f()
      f_calc = fc_sq.f_sq_as_f()
      R1 = f_obs.r1_factor(f_calc,
        scale_factor=math.sqrt(self.scale_factor()), assume_index_matching=True)
      return R1, f_obs.size()

    def covariance_matrix(self,
                          jacobian_transpose=None,
                          normalised_by_goof=True):
      """ The columns of the jacobian_transpose determine which crystallographic
          parameters appear in the covariance matrix.
          If jacobian_transpose is None, then the covariance matrix returned will
          be that for the independent L.S. parameters.
      """
      if not getattr(self, 'normal_equations_are_complete', True):
        """ Build what is missing rather than refusing.

        Conjugate gradients does not form the normal equations, and a run
        which declined the standard uncertainties closed on an objective-only
        pass -- so what is standing holds the restraints and the origin fixing
        and nothing from the data. Inverting that fails in Cholesky on the
        first parameter no restraint happens to touch, which says nothing
        about that parameter and has sent more than one investigation after an
        innocent atom.

        Declining the s.u. is meant to save the closing build, not to make the
        covariance matrix unobtainable. So it is built here, at the parameters
        the refinement finished on -- which is exactly what the closing build
        would have done had the s.u. been asked for in the first place. One
        pass over the reflections, paid only by a caller that wants a
        covariance matrix, and paid once because the flag is then set.
        """
        # weighted with the scale the solver finished on, not the one the
        # opening pass left behind; build_up clears the override itself
        solver_scale = getattr(self, 'solver_scale_factor', None)
        if solver_scale is not None and self.overridden_scale_factor is None:
          self.overridden_scale_factor = solver_scale
        self.build_up()
        if not self.normal_equations_are_complete:
          raise RuntimeError(
            "the normal equations are still incomplete after a full build; "
            "there is no covariance matrix to take.")
      if not self.step_equations().solved:
        self.solve()
      cov = linalg.inverse_of_u_transpose_u(
        self.step_equations().cholesky_factor_packed_u())
      cov /= self.sum_w_yo_sq()
      if jacobian_transpose is not None:
        cov = jacobian_transpose.self_transpose_times_symmetric_times_self(cov)
      if normalised_by_goof: cov *= self.variance_goof_factor()
      return cov

    def covariance_matrix_and_annotations(self):
      jac_tr = self.reparametrisation.jacobian_transpose_matching_grad_fc()
      return covariance_matrix_and_annotations(
        self.covariance_matrix(jacobian_transpose=jac_tr),
        self.reparametrisation.component_annotations)

  return klass


def crystallographic_ls(
  observations, reparametrisation,
  non_linear_ls_with_separable_scale_factor=None,
  may_parallelise=True,
  **kwds):
  return crystallographic_ls_class(
    non_linear_ls_with_separable_scale_factor,
    n_parameters=reparametrisation.n_independents,
    may_parallelise=may_parallelise)(observations, reparametrisation,
                                     may_parallelise=may_parallelise, **kwds)


class covariance_matrix_and_annotations(object):

  def __init__(self, covariance_matrix, annotations):
    """ The covariance matrix is assumed to be a symmetric matrix stored as a
        packed upper diagonal matrix.
    """
    self.matrix = covariance_matrix
    self.annotations = annotations
    self._2_n_minus_1 = 2*len(self.annotations)-1 # precompute for efficiency

  def __call__(self, i, j):
    return self.matrix[i*(self._2_n_minus_1-i)//2 + j]

  def variance_of(self, annotation):
    i = self.annotations.index(annotation)
    return self(i, i)

  def covariance_of(self, annotation_1, annotation_2):
    i = self.annotations.index(annotation_1)
    j = self.annotations.index(annotation_2)
    if j > i:
      i, j = j, i
    return self(i, j)

  def diagonal(self):
    return self.matrix.matrix_packed_u_diagonal()
