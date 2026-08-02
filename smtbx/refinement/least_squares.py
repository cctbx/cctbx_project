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
                              may_parallelise=None):
  """ Construct a class for crystallographic L.S. based on the given engine

  n_parameters and may_parallelise are the problem this class is about to be
  used on, where the caller knows them. They only ever pick between the two
  stock accumulators; an engine passed in explicitly is always honoured, and
  with neither given the choice is exactly what it has always been.
  """
  def get_base_class(non_linear_ls_with_separable_scale_factor):
    base_class = non_linear_ls_with_separable_scale_factor
    if not base_class:
      try:
        from fast_linalg import env
        if env.initialised:
          base_class = normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_3
        else:
          base_class = normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2
      except Exception:
        base_class = normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2
      # A threaded build of a small enough problem is the one case where the
      # packed rank-1 update beats the buffered rank-k one. Each thread
      # accumulates into a matrix of its own, and while that matrix is small it
      # stays in cache, where BLAS 3's row buffer, its n x n scratch and the
      # per-thread syrks competing with each other all cost more than they save.
      if (may_parallelise
          and n_parameters is not None
          and n_parameters <= blas_2_parallel_max_parameters()):
        base_class = normal_eqns.non_linear_ls_with_separable_scale_factor_BLAS_2
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

      def build_normal_eqns(scale_factor, weighting_scheme, objective_only):
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
          self.max_memory)

      if not self.finalised: #i.e. never been called
        self.reparametrisation.linearise()
        self.reparametrisation.store()
        scale_factor = self.initial_scale_factor
        if scale_factor is None: # we haven't got one from previous refinement
          result = build_normal_eqns(scale_factor=None,
                                     weighting_scheme=sigma_weighting(),
                                     objective_only=True)
          scale_factor = self.scale_factor()
      else: # use scale factor from the previous step
        scale_factor = self.scale_factor()

      self.reset()
      result = build_normal_eqns(scale_factor,
                                 self.weighting_scheme,
                                 objective_only)
      self.f_calc = self.observations.fo_sq.array(
        data=result.f_calc(), sigmas=None)
      self.fc_sq = self.observations.fo_sq.array(
        data=result.observables(), sigmas=None)\
          .set_observation_type_xray_intensity()
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
      # An override says which scale factor to weight the data with when there
      # are no fresh normal equations to take one from. There now are, so it
      # has served its purpose, and anything reading scale_factor() after this
      # -- the journal in normal_eqns_solving, for one -- wants the new one.
      self.overridden_scale_factor = None

    def parameter_vector_norm(self):
      return self.reparametrisation.norm_of_independent_parameter_vector

    # Set by a solver which determines the scale factor without building the
    # normal equations, c.f. smtbx.refinement.cgls; None means take it from
    # them as usual. build_up clears it, having made it obsolete.
    overridden_scale_factor = None

    def scale_factor(self):
      if self.overridden_scale_factor is not None:
        return self.overridden_scale_factor
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
      return math.sqrt(self.chi_sq_data_only)

    def restrained_goof(self):
      if self.restraints_manager is None:
        return self.goof()
      return math.sqrt(self.chi_sq_data_and_restraints)

    def wR2(self, cutoff_factor=None):
      if cutoff_factor is None:
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
      if not self.step_equations().solved:
        self.solve()
      cov = linalg.inverse_of_u_transpose_u(
        self.step_equations().cholesky_factor_packed_u())
      cov /= self.sum_w_yo_sq()
      if jacobian_transpose is not None:
        cov = jacobian_transpose.self_transpose_times_symmetric_times_self(cov)
      if normalised_by_goof: cov *= self.restrained_goof()**2
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
