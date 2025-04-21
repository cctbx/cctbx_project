from __future__ import absolute_import, division, print_function
import six
from six.moves import range
from six.moves import cStringIO as StringIO
from collections import Counter
import math
from xfel.merging.application.worker import worker
from libtbx import adopt_init_args, group_args
from dials.array_family import flex
from dxtbx.model.experiment_list import ExperimentList
from cctbx import miller
from cctbx.crystal import symmetry
from scitbx import matrix
from scitbx.math.tests.tst_weighted_correlation import simple_weighted_correlation
from cctbx.crystal_orientation import crystal_orientation, basis_type

class postrefinement_rs(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(postrefinement_rs, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

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
      I_observed = observations_pair1_selected.data()
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

      if self.params.postrefinement.algorithm == "rs":

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
        refinery = rs_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE, ICALCVEC = I_reference, IOBSVEC = I_observed)

      elif self.params.postrefinement.algorithm == "eta_deff":

        eta_init = 2. * MOSAICITY_DEG * math.pi/180.
        D_eff_init = 2. * DOMAIN_SIZE_A
        current = flex.double([SWC.slope, BFACTOR, eta_init, 0., 0., D_eff_init])

        parameterization_class = eta_deff_parameterization
        refinery = eta_deff_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE, ICALCVEC = I_reference, IOBSVEC = I_observed)

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
        dcl = self.params.postrefinement.delta_corr_limit
        if dcl is not None and SWC_after_post.corr - SWC.corr < -(dcl):
          raise ValueError("CC decreased by > {} in postrefinement".format(dcl))
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
          if not key in new_exp_reflections.keys() and key in exp_reflections_match_results.keys():
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
    self.MINI = lbfgs_minimizer_base(self.params,
                                     current_x = self.current,
                                     parameterization = self.parameterization_class,
                                     refinery = self.refinery,
                                     out = out)
    if self.params.output.log_level == 0:
      self.logger.log("\n" + out.getvalue())

  def result_for_cxi_merge(self):

    scaler = self.refinery.scaler_callable(self.parameterization_class(self.MINI.x))

    if self.params.postrefinement.algorithm == "rs":
      fat_selection = (self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x)) > self.params.postrefinement.partiality_threshold_hcfix)
      fats = self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x))
    else:
      fat_selection = (self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x)) < 0.9)

    fat_count = fat_selection.count(True)

    # reject an experiment with insufficient number of near-full reflections
    if fat_count < 3:

      if self.params.output.log_level == 0:
        self.logger.log("Rejected experiment, because: On total %5d the fat selection is %5d"%(len(self.observations_pair1_selected.indices()), fat_count))

      '''
      # debugging
      rejected_fat_max = 0.0
      for fat in fats:
        if fat <= 0.2:
          if fat > rejected_fat_max:
            rejected_fat_max = fat
      self.logger.log("MAXIMUM FAT VALUE AMONG REJECTED REFLECTIONS IS: %f"%rejected_fat_max)
      '''

      raise ValueError("< 3 near-fulls after refinement")

    if self.params.output.log_level == 0:
      self.logger.log("On total %5d the fat selection is %5d"%(len(self.observations_pair1_selected.indices()), fat_count))

    observations_original_index = self.observations_original_index_pair1_selected.select(fat_selection)

    observations = self.observations_pair1_selected.customized_copy(indices = self.observations_pair1_selected.indices().select(fat_selection),
                                                                    data = (self.observations_pair1_selected.data()/scaler).select(fat_selection),
                                                                    sigmas = (self.observations_pair1_selected.sigmas()/scaler).select(fat_selection))

    matches = miller.match_multi_indices(miller_indices_unique = self.params.scaling.miller_set.indices(),
                                         miller_indices = observations.indices())

    return observations_original_index, observations, matches

  def get_parameter_values(self):
    values = self.parameterization_class(self.MINI.x)
    return values

  def result_for_samosa(self):
    values = self.parameterization_class(self.MINI.x)
    return self.refinery.get_eff_Astar(values), values.RS

class refinery_base(group_args):
    def __init__(self, **kwargs):
      group_args.__init__(self,**kwargs)
      mandatory = ["ORI","MILLER","BEAM","WAVE","ICALCVEC","IOBSVEC"]
      for key in mandatory: getattr(self,key)
      self.DSSQ = self.ORI.unit_cell().d_star_sq(self.MILLER)

    """Refinery class takes reference and observations, and implements target
    functions and derivatives for a particular model paradigm."""
    def get_Rh_array(self, values):
      eff_Astar = self.get_eff_Astar(values)
      h = self.MILLER.as_vec3_double()
      x = flex.mat3_double(len(self.MILLER), eff_Astar) * h
      Svec = x + self.BEAM
      Rh = Svec.norms() - (1./self.WAVE)
      return Rh

    def get_s1_array(self, values):
      miller_vec = self.MILLER.as_vec3_double()
      ref_ori = matrix.sqr(self.ORI.reciprocal_matrix())
      Rx = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(values.thetax)
      Ry = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(values.thetay)
      s_array = flex.mat3_double(len(self.MILLER),Ry * Rx * ref_ori) * miller_vec
      s1_array = s_array + flex.vec3_double(len(self.MILLER), self.BEAM)
      return s1_array

    def get_eff_Astar(self, values):
      thetax = values.thetax; thetay = values.thetay;
      effective_orientation = self.ORI.rotate_thru((1,0,0),thetax
         ).rotate_thru((0,1,0),thetay
         )
      return matrix.sqr(effective_orientation.reciprocal_matrix())

    def scaler_callable(self, values):
      PB = self.get_partiality_array(values)
      EXP = flex.exp(-2.*values.BFACTOR*self.DSSQ)
      terms = values.G * EXP * PB
      return terms

    def fvec_callable(self, values):
      PB = self.get_partiality_array(values)
      EXP = flex.exp(-2.*values.BFACTOR*self.DSSQ)
      terms = (values.G * EXP * PB * self.ICALCVEC - self.IOBSVEC)
      # Ideas for improvement
      #   straightforward to also include sigma weighting
      #   add extra terms representing rotational excursion: terms.concatenate(1.e7*Rh)
      return terms

class rs_refinery(refinery_base):
    def lorentz_callable(self,values):
      return self.get_partiality_array(values)

    def get_partiality_array(self,values):
      rs = values.RS
      Rh = self.get_Rh_array(values)
      rs_sq = rs*rs
      PB = rs_sq / ((2. * (Rh * Rh)) + rs_sq)
      return PB

class eta_deff_refinery(refinery_base):
    def __init__(self, **kwargs):
      refinery_base.__init__(self,**kwargs)
      self.DVEC = self.ORI.unit_cell().d(self.MILLER)

    def lorentz_callable(self,values):
      Rh = self.get_Rh_array(values)
      Rs = flex.double(len(self.MILLER),1./values.DEFF)+flex.double(len(self.MILLER),values.ETA/2.)/self.DVEC
      ratio = Rh / Rs
      ratio_abs = flex.abs(ratio)
      return ratio_abs

    def get_partiality_array(self,values):
      Rh = self.get_Rh_array(values)
      Rs = flex.double(len(self.MILLER),1./values.DEFF)+flex.double(len(self.MILLER),values.ETA/2.)/self.DVEC
      Rs_sq = Rs * Rs
      Rh_sq = Rh * Rh
      numerator = Rs_sq - Rh_sq
      denominator = values.DEFF * Rs * Rs_sq
      partiality = numerator / denominator
      return partiality

class unpack_base(object):
  "abstract interface"
  def __init__(YY,values):
    YY.reference = values # simply the flex double list of parameters
    YY.keys = None # override
  def __getattr__(YY,item):
    if item not in YY.keys:
      raise AttributeError(item)
    return YY.reference[YY.keys.index(item)]
  def show(values,out):
    raise NotImplementedError

class rs_parameterization(unpack_base):
  def __init__(YY,values):
    super(rs_parameterization, YY).__init__(values)
    YY.keys = ['G', 'BFACTOR','RS','thetax','thetay']

  def show(YY, out):
    print ("G: %10.7f; B: %10.7f; RS: %10.7f; THETAX: %7.3f deg; THETAY: %7.3f deg"\
          %(YY.G, YY.BFACTOR, YY.RS, 180.*YY.thetax/math.pi, 180.*YY.thetay/math.pi)\
          , file=out)

class eta_deff_parameterization(unpack_base):
  def __init__(YY,values):
    super(eta_deff_parameterization, YY).__init__(values)
    YY.keys = ['G', 'BFACTOR','ETA','thetax','thetay','DEFF']

  def show(YY, out):
    print ("G: %10.7f; B: %10.7f; eta: %10.7f; Deff %10.2f; THETAX: %7.3f deg; THETAY: %7.3f deg"\
          %(YY.G, YY.BFACTOR, YY.ETA, YY.DEFF, 180.*YY.thetax/math.pi, 180.*YY.thetay/math.pi)\
          , file=out)

class lbfgs_minimizer_base:

  def __init__(self, params, current_x=None, parameterization=None, refinery=None, out=None,
               min_iterations=0, max_calls=1000, max_drop_eps=1.e-5):
    adopt_init_args(self, locals())
    self.n = current_x.size()
    self.x = current_x
    from scitbx import lbfgs
    self.minimizer = lbfgs.run(
      target_evaluator=self,
      termination_params=lbfgs.termination_parameters(
        traditional_convergence_test=False,
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
    assert -150. < values.BFACTOR < 150. # limits on the exponent, please
    self.func = self.refinery.fvec_callable(values)
    functional = flex.sum(self.func*self.func)
    self.f = functional
    DELTA = 1.E-7
    self.g = flex.double()
    for x in range(self.n):
      templist = list(self.x)
      templist[x]+=DELTA
      dvalues = flex.double(templist)

      dfunc = self.refinery.fvec_callable(self.parameterization(dvalues))
      dfunctional = flex.sum(dfunc*dfunc)
      #calculate by finite_difference
      self.g.append( ( dfunctional-functional )/DELTA )

    if self.params.postrefinement.algorithm == 'rs':
      for p in self.params.postrefinement.rs.fix:
        self.g[values.keys.index(p)] = 0

    print ("rms %10.3f; "%math.sqrt(flex.mean(self.func*self.func)), file=self.out, end='')
    values.show(self.out)

    return self.f, self.g

  def __del__(self):
    values = self.parameterization(self.x)
    print ("FINALMODEL", file=self.out)
    print ("rms %10.3f; "%math.sqrt(flex.mean(self.func*self.func)), file=self.out, end='')
    values.show(self.out)

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(postrefinement_rs)
