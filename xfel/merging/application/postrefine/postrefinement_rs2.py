from __future__ import absolute_import, division, print_function
import six
from six.moves import range
from six.moves import cStringIO as StringIO
from collections import Counter
import math
from libtbx import adopt_init_args
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from cctbx import miller
from cctbx.crystal import symmetry
from scitbx import matrix
from scitbx.math.tests.tst_weighted_correlation import simple_weighted_correlation
from cctbx.crystal_orientation import crystal_orientation, basis_type
from xfel.merging.application.postrefine.postrefinement_rs import postrefinement_rs, rs_refinery, rs_parameterization, lbfgs_minimizer_base

def chosen_weights(observation_set, params):
    data = observation_set.data()
    sigmas = observation_set.sigmas()
    return {
      "unit": flex.double(len(data),1.),
      "variance": 1./(sigmas*sigmas),
      "gentle": flex.pow(flex.sqrt(flex.abs(data))/sigmas,2),
      "extreme": flex.pow(data/sigmas,2)
    } [ params.postrefinement.target_weighting ]

class postrefinement_rs2(postrefinement_rs):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(postrefinement_rs2, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Postrefinement'

  def run(self, experiments, reflections):
    self.logger.log_step_time("POSTREFINEMENT")
    if (not self.params.postrefinement.enable) or (self.params.scaling.algorithm != "mark0"): # mark1 implies no scaling/post-refinement
      self.logger.log("No post-refinement was done")
      if self.mpi_helper.rank == 0:
        self.logger.main_log("No post-refinement was done")
      return experiments, reflections

    target_symm = symmetry(unit_cell = self.params.scaling.unit_cell, space_group_info = self.params.scaling.space_group)
    i_model = self.params.scaling.i_model
    miller_set = self.params.scaling.miller_set

    # Ensure that match_multi_indices() will return identical results
    # when a frame's observations are matched against the
    # pre-generated Miller set, miller_set, and the reference
    # data set, i_model.  The implication is that the same match
    # can be used to map Miller indices to array indices for intensity
    # accumulation, and for determination of the correlation
    # coefficient in the presence of a scaling reference.
    assert len(i_model.indices()) == len(miller_set.indices())
    assert (i_model.indices() == miller_set.indices()).count(False) == 0

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()

    experiments_rejected_by_reason = Counter()  # reason:how_many_rejected

    for expt_id, experiment in enumerate(experiments):

      exp_reflections = reflections.select(reflections['id'] == expt_id)

      # Build a miller array with _original_ miller indices of the experiment reflections
      exp_miller_indices_original = miller.set(target_symm, exp_reflections['miller_index'], not self.params.merging.merge_anomalous)
      observations_original_index = miller.array(exp_miller_indices_original, exp_reflections['intensity.sum.value'], flex.sqrt(exp_reflections['intensity.sum.variance']))

      assert exp_reflections.size() == exp_miller_indices_original.size()
      assert observations_original_index.size() == exp_miller_indices_original.size()

      # Build a miller array with _asymmetric_ miller indices of the experiment reflections
      exp_miller_indices_asu = miller.set(target_symm, exp_reflections['miller_index_asymmetric'], True)
      observations = miller.array(exp_miller_indices_asu, exp_reflections['intensity.sum.value'], flex.sqrt(exp_reflections['intensity.sum.variance']))

      matches = miller.match_multi_indices(miller_indices_unique = miller_set.indices(), miller_indices = observations.indices())

      pair1 = flex.int([pair[1] for pair in matches.pairs()]) # refers to the observations
      pair0 = flex.int([pair[0] for pair in matches.pairs()]) # refers to the model

      # narrow things down to the set that matches, only
      observations_pair1_selected = observations.customized_copy(indices = flex.miller_index([observations.indices()[p] for p in pair1]),
                                                                 data = flex.double([observations.data()[p] for p in pair1]),
                                                                 sigmas = flex.double([observations.sigmas()[p] for p in pair1]))

      observations_original_index_pair1_selected = observations_original_index.customized_copy(indices = flex.miller_index([observations_original_index.indices()[p] for p in pair1]),
                                                                                               data = flex.double([observations_original_index.data()[p] for p in pair1]),
                                                                                               sigmas = flex.double([observations_original_index.sigmas()[p] for p in pair1]))
###################
      I_observed = observations_pair1_selected.data()
      chosen = chosen_weights(observations_pair1_selected, self.params)

      MILLER = observations_original_index_pair1_selected.indices()
      ORI = crystal_orientation(experiment.crystal.get_A(), basis_type.reciprocal)
      Astar = matrix.sqr(ORI.reciprocal_matrix())
      Astar_from_experiment = matrix.sqr(experiment.crystal.get_A())
      assert Astar == Astar_from_experiment

      WAVE = experiment.beam.get_wavelength()
      BEAM = matrix.col((0.0,0.0,-1./WAVE))
      BFACTOR = 0.
      MOSAICITY_DEG = experiment.crystal.get_half_mosaicity_deg()
      DOMAIN_SIZE_A = experiment.crystal.get_domain_size_ang()

      # calculation of correlation here
      I_reference = flex.double([i_model.data()[pair[0]] for pair in matches.pairs()])
      I_invalid = flex.bool([i_model.sigmas()[pair[0]] < 0. for pair in matches.pairs()])
      use_weights = False # New facility for getting variance-weighted correlation

      if use_weights:
        # variance weighting
        I_weight = flex.double([1./(observations_pair1_selected.sigmas()[pair[1]])**2 for pair in matches.pairs()])
      else:
        I_weight = flex.double(len(observations_pair1_selected.sigmas()), 1.)

      I_weight.set_selected(I_invalid, 0.)
      chosen.set_selected(I_invalid,0.)

      """Explanation of 'include_negatives' semantics as originally implemented in cxi.merge postrefinement:
         include_negatives = True
         + and - reflections both used for Rh distribution for initial estimate of RS parameter
         + and - reflections both used for calc/obs correlation slope for initial estimate of G parameter
         + and - reflections both passed to the refinery and used in the target function (makes sense if
                             you look at it from a certain point of view)

         include_negatives = False
         + and - reflections both used for Rh distribution for initial estimate of RS parameter
         +       reflections only used for calc/obs correlation slope for initial estimate of G parameter
         + and - reflections both passed to the refinery and used in the target function (makes sense if
                             you look at it from a certain point of view)

         NOTE: by the new design, "include negatives" is always True
      """

      SWC = simple_weighted_correlation(I_weight, I_reference, I_observed)
      if self.params.output.log_level == 0:
        self.logger.log("Old correlation is: %f"%SWC.corr)

      Rhall = flex.double()
      for mill in MILLER:
        H = matrix.col(mill)
        Xhkl = Astar*H
        Rh = ( Xhkl + BEAM ).length() - (1./WAVE)
        Rhall.append(Rh)
      Rs = math.sqrt(flex.mean(Rhall*Rhall))

      RS = 1./10000. # reciprocal effective domain size of 1 micron
      RS = Rs        # try this empirically determined approximate, monochrome, a-mosaic value
      current = flex.double([SWC.slope, BFACTOR, RS, 0., 0.])

      parameterization_class = rs_parameterization
      refinery = rs2_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE, ICALCVEC = I_reference, IOBSVEC = I_observed, WEIGHTS = chosen)
      refinery.set_profile_shape(self.params.postrefinement.lineshape)

      func = refinery.fvec_callable(parameterization_class(current))
      functional = flex.sum(func * func)

      if self.params.output.log_level == 0:
        self.logger.log("functional: %f"%functional)

      self.current = current;
      self.parameterization_class = parameterization_class
      self.refinery = refinery;

      self.observations_pair1_selected = observations_pair1_selected;
      self.observations_original_index_pair1_selected = observations_original_index_pair1_selected

      error_detected = False

      try:
        self.run_plain()

        result_observations_original_index, result_observations, result_matches = self.result_for_cxi_merge()

        assert result_observations_original_index.size() == result_observations.size()
        assert result_matches.pairs().size() == result_observations_original_index.size()
        # Calculate the correlation of each frame after corrections.
        # This is used in the MM24 error model to determine a per frame level of error
        # These are the added to the reflection table
        I_observed = result_observations.data()
        matches = miller.match_multi_indices(
          miller_indices_unique = miller_set.indices(),
          miller_indices = result_observations.indices()
        )
        I_reference = flex.double([i_model.data()[pair[0]] for pair in matches.pairs()])
        I_invalid = flex.bool([i_model.sigmas()[pair[0]] < 0. for pair in matches.pairs()])
        I_weight = flex.double(len(result_observations.sigmas()), 1.)
        I_weight.set_selected(I_invalid, 0.)
        SWC_after_post = simple_weighted_correlation(I_weight, I_reference, I_observed)
      except (AssertionError, ValueError, RuntimeError) as e:
        error_detected = True
        reason = repr(e)
        if not reason:
          reason = "Unknown error"
        experiments_rejected_by_reason[reason] += 1

      if not error_detected:
        new_experiments.append(experiment)

        new_exp_reflections = flex.reflection_table()
        new_exp_reflections['miller_index_asymmetric']  = result_observations.indices()
        new_exp_reflections['intensity.sum.value']      = result_observations.data()
        new_exp_reflections['intensity.sum.variance']   = flex.pow(result_observations.sigmas(),2)
        new_exp_reflections['id']                       = flex.int(len(new_exp_reflections), len(new_experiments)-1)
        new_exp_reflections.experiment_identifiers()[len(new_experiments)-1] = experiment.identifier

        # The original reflection table, i.e. the input to this run() method, has more columns than those used
        # for the postrefinement ("data" and "sigma" in the miller arrays). The problems is: some of the input reflections may have been rejected by now.
        # So to bring those extra columns over to the new reflection table, we have to create a subset of the original exp_reflections table,
        # which would match (by original miller indices) the miller array results of the postrefinement.
        match_original_indices = miller.match_multi_indices(miller_indices_unique = exp_miller_indices_original.indices(), miller_indices = result_observations_original_index.indices())
        exp_reflections_match_results = exp_reflections.select(match_original_indices.pairs().column(0))
        assert (exp_reflections_match_results['intensity.sum.value'] == result_observations_original_index.data()).count(False) == 0
        new_exp_reflections['intensity.sum.value.unmodified'] = exp_reflections_match_results['intensity.sum.value.unmodified']
        new_exp_reflections['intensity.sum.variance.unmodified'] = exp_reflections_match_results['intensity.sum.variance.unmodified']
        for key in self.params.input.persistent_refl_cols:
          if key not in new_exp_reflections.keys():
            new_exp_reflections[key] = exp_reflections_match_results[key]
        new_exp_reflections["correlation_after_post"] = flex.double(len(new_exp_reflections), SWC_after_post.corr)
        new_reflections.extend(new_exp_reflections)

    # report rejected experiments, reflections
    experiments_rejected_by_postrefinement = len(experiments) - len(new_experiments)
    reflections_rejected_by_postrefinement = reflections.size() - new_reflections.size()

    self.logger.log("Experiments rejected by post-refinement: %d"%experiments_rejected_by_postrefinement)
    self.logger.log("Reflections rejected by post-refinement: %d"%reflections_rejected_by_postrefinement)

    for reason, count in six.iteritems(experiments_rejected_by_reason):
      self.logger.log("Experiments rejected due to %s: %d"%(reason,count))

    # Now that each rank has all reasons from all ranks, we can treat the reasons in a uniform way.
    total_experiments_rejected_by_reason = self.mpi_helper.count(experiments_rejected_by_reason)
    total_accepted_experiment_count = self.mpi_helper.sum(len(new_experiments))

    # how many reflections have we rejected due to post-refinement?
    rejected_reflections = len(reflections) - len(new_reflections);
    total_rejected_reflections = self.mpi_helper.sum(rejected_reflections)

    if self.mpi_helper.rank == 0:
      for reason, count in six.iteritems(total_experiments_rejected_by_reason):
        self.logger.main_log("Total experiments rejected due to %s: %d"%(reason,count))
      self.logger.main_log("Total experiments accepted: %d"%total_accepted_experiment_count)
      self.logger.main_log("Total reflections rejected due to post-refinement: %d"%total_rejected_reflections)

    self.logger.log_step_time("POSTREFINEMENT", True)

    # Do we have any data left?
    from xfel.merging.application.utils.data_counter import data_counter
    data_counter(self.params).count(new_experiments, new_reflections)

    return new_experiments, new_reflections

  def run_plain(self):

    out = StringIO()
    self.MINI = lbfgs_minimizer_derivatives(self.params,
                                           current_x = self.current,
                                           parameterization = self.parameterization_class,
                                           refinery = self.refinery,
                                           out = out)
    self.refined_mini = self.MINI

    if self.params.output.log_level == 0:
      self.logger.log("\n" + out.getvalue())

  def rs2_parameter_range_assertions(self, values):
    # New range assertions for refined variables
    assert 0 < values.G, "G-scale value out of range ( < 0 ) after rs2 refinement"
    assert -25 < values.BFACTOR and values.BFACTOR < 25, "B-factor value out of range ( |B|>25 ) after rs2 refinement"
    assert -0.5<180.*values.thetax/math.pi<0.5,"thetax value out of range ( |rotx|>.5 degrees ) after rs2 refinement"
    assert -0.5<180.*values.thetay/math.pi<0.5,"thetay value out of range ( |roty|>.5 degrees ) after rs2 refinement"

  def result_for_cxi_merge(self):
    values = self.get_parameter_values()
    self.rs2_parameter_range_assertions(values)
    scaler = self.refinery.scaler_callable(self.parameterization_class(self.MINI.x))

    partiality_array = self.refinery.get_partiality_array(values)
    p_scaler = flex.pow(partiality_array,
                        0.5*self.params.postrefinement.merge_partiality_exponent)
    fat_selection = (partiality_array > self.params.postrefinement.partiality_threshold_hcfix)
    fat_count = fat_selection.count(True)
    scaler_s = scaler.select(fat_selection)
    p_scaler_s = p_scaler.select(fat_selection)

    # reject an experiment with insufficient number of near-full reflections
    if fat_count < 3:
      if self.params.output.log_level == 0:
        self.logger.log("Rejected experiment, because: On total %5d the fat selection is %5d"%(len(self.observations_pair1_selected.indices()), fat_count))
      raise ValueError("< 3 near-fulls after refinement")
    if self.params.output.log_level == 0:
      self.logger.log("On total %5d the fat selection is %5d"%(len(self.observations_pair1_selected.indices()), fat_count))

    observations_original_index = self.observations_original_index_pair1_selected.select(fat_selection)

    observations = self.observations_pair1_selected.customized_copy(
      indices = self.observations_pair1_selected.indices().select(fat_selection),
      data = (self.observations_pair1_selected.data().select(fat_selection)/scaler_s),
      sigmas = (self.observations_pair1_selected.sigmas().select(fat_selection)/(scaler_s * p_scaler_s))
    )
    matches = miller.match_multi_indices(
      miller_indices_unique=self.params.scaling.miller_set.indices(),
      miller_indices=observations.indices())

    I_weight = flex.double(len(observations.sigmas()), 1.)
    I_reference = flex.double([self.params.scaling.i_model.data()[pair[0]] for pair in matches.pairs()])
    I_invalid = flex.bool([self.params.scaling.i_model.sigmas()[pair[0]] < 0. for pair in matches.pairs()])
    I_weight.set_selected(I_invalid,0.)
    SWC = simple_weighted_correlation(I_weight, I_reference, observations.data())

    if self.params.output.log_level == 0:
      self.logger.log("CORR: NEW correlation is: %f"%SWC.corr)
      self.logger.log("ASTAR: ")
      self.logger.log(tuple(self.refinery.get_eff_Astar(values)))

    self.final_corr = SWC.corr
    self.refined_mini = self.MINI

    #another range assertion
    assert self.final_corr > 0.1,"correlation coefficient out of range (<= 0.1) after rs2 refinement"

    return observations_original_index, observations, matches

  def get_parameter_values(self):
    return self.refined_mini.parameterization(self.refined_mini.x)

class rs2_refinery(rs_refinery):

    def set_profile_shape(self, shape):
      self.profile_shape = shape
      self.get_partiality_array = {
        "lorentzian":super(rs2_refinery, self).get_partiality_array,
        "gaussian": self.get_gaussian_partiality_array
      }[shape]

    def get_gaussian_partiality_array(self,values):
      rs = values.RS
      Rh = self.get_Rh_array(values)
      immersion = Rh/rs
      gaussian = flex.exp(-2. * math.log(2) * (immersion*immersion))
      return gaussian

    def jacobian_callable(self,values):
      PB = self.get_partiality_array(values)
      EXP = flex.exp(-2.*values.BFACTOR*self.DSSQ)
      G_terms = (EXP * PB * self.ICALCVEC)
      B_terms = (values.G * EXP * PB * self.ICALCVEC)*(-2.*self.DSSQ)
      P_terms = (values.G * EXP * self.ICALCVEC)

      thetax = values.thetax; thetay = values.thetay;
      Rx = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(thetax)
      dRx_dthetax = matrix.col((1,0,0)).axis_and_angle_as_r3_derivative_wrt_angle(thetax)
      Ry = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(thetay)
      dRy_dthetay = matrix.col((0,1,0)).axis_and_angle_as_r3_derivative_wrt_angle(thetay)
      ref_ori = matrix.sqr(self.ORI.reciprocal_matrix())
      miller_vec = self.MILLER.as_vec3_double()
      ds1_dthetax = flex.mat3_double(len(self.MILLER),Ry * dRx_dthetax * ref_ori) * miller_vec
      ds1_dthetay = flex.mat3_double(len(self.MILLER),dRy_dthetay * Rx * ref_ori) * miller_vec

      s1vec = self.get_s1_array(values)
      s1lenvec = flex.sqrt(s1vec.dot(s1vec))
      dRh_dthetax = s1vec.dot(ds1_dthetax)/s1lenvec
      dRh_dthetay = s1vec.dot(ds1_dthetay)/s1lenvec
      rs = values.RS
      Rh = self.get_Rh_array(values)
      rs_sq = rs*rs
      denomin = (2. * Rh * Rh + rs_sq)
      dPB_dRh = { "lorentzian": -PB * 4. * Rh / denomin,
                  "gaussian": -PB * 4. * math.log(2) * Rh / rs_sq }[self.profile_shape]
      dPB_dthetax = dPB_dRh * dRh_dthetax
      dPB_dthetay = dPB_dRh * dRh_dthetay
      Px_terms = P_terms * dPB_dthetax; Py_terms = P_terms * dPB_dthetay

      return [G_terms,B_terms,0,Px_terms,Py_terms]

class lbfgs_minimizer_derivatives(lbfgs_minimizer_base):

  def __init__(self, params, current_x=None, parameterization=None, refinery=None, out=None,
               min_iterations=0, max_calls=1000, max_drop_eps=1.e-5):
    adopt_init_args(self, locals())
    self.n = current_x.size()
    self.x = current_x
    from scitbx import lbfgs
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs.termination_parameters(
        traditional_convergence_test=True,
        drop_convergence_test_max_drop_eps=max_drop_eps,
        min_iterations=min_iterations,
        max_iterations = None,
        max_calls=max_calls),
      exception_handling_params=lbfgs.exception_handling_parameters(
         ignore_line_search_failed_rounding_errors=True,
         ignore_line_search_failed_step_at_lower_bound=True,#the only change from default
         ignore_line_search_failed_step_at_upper_bound=False,
         ignore_line_search_failed_maxfev=False,
         ignore_line_search_failed_xtol=False,
         ignore_search_direction_not_descent=False)
      )

  def compute_functional_and_gradients(self):
    values = self.parameterization(self.x)
    assert -150. < values.BFACTOR < 150,"B-factor out of range (+/-150) within rs2 functional and gradients"
    self.func = self.refinery.fvec_callable(values)
    functional = flex.sum(self.refinery.WEIGHTS*self.func*self.func)
    self.f = functional
    jacobian = self.refinery.jacobian_callable(values)
    self.g = flex.double(self.n)
    for ix in range(self.n):
      self.g[ix] = flex.sum(2. * self.refinery.WEIGHTS * self.func * jacobian[ix])

    print("rms %10.3f"%math.sqrt(flex.sum(self.refinery.WEIGHTS*self.func*self.func)/flex.sum(self.refinery.WEIGHTS)), file=self.out, end='')

    values.show(self.out)

    return self.f, self.g

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(postrefinement_rs2)
