from __future__ import absolute_import, division, print_function
from six.moves import range
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
from six.moves import cStringIO as StringIO
from libtbx.development.timers import Profiler
from libtbx.test_utils import approx_equal
import six

class ratchet1_postrefinement(worker):

  def __init__(self, params, mpi_helper=None, mpi_logger=None):
    super(ratchet1_postrefinement, self).__init__(params=params, mpi_helper=mpi_helper, mpi_logger=mpi_logger)

  def __repr__(self):
    return 'Ratchet1_postrefinement'

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
    # pre-generated Miller set, self.miller_set, and the reference
    # data set, self.i_model.  The implication is that the same match
    # can be used to map Miller indices to array indices for intensity
    # accumulation, and for determination of the correlation
    # coefficient in the presence of a scaling reference.
    assert len(i_model.indices()) == len(miller_set.indices())
    assert (i_model.indices() == miller_set.indices()).count(False) == 0

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()
    # seems like we'll need a backwards mapping from refl to expt.  Here it is:
    # experiments.find(reflections[refl_id]['exp_id'])

    experiments_rejected_by_reason = {} # reason:how_many_rejected

    # deep dive to create new reference columns in the input table
    # doesn't look like it will work if asu indices are missing in the i_model
    # flex.<miller_index>.set_selected() is a very powerful tool
    # M.set_selected(S,U): S a flex.bool of same length M, U a single value
    # M.set_selected(S,U): S a flex.bool of same length M, U a flex array of same length as # of True values
    # M.set_selected(S,U): S a flex.size_t of addresses into M, U a flex array of same length as S

    new_matches = miller.match_multi_indices(miller_indices_unique = miller_set.indices(),
                                             miller_indices = reflections['miller_index_asymmetric'])
    pair1 = flex.size_t([pair[1] for pair in new_matches.pairs()]) # refers to the observations
    pair0 = flex.size_t([pair[0] for pair in new_matches.pairs()]) # refers to the model
    IR = intensity_reference = flex.double(len(reflections))

    matched_intensity_references = self.params.scaling.i_model.data().select(pair0)
    intensity_reference.set_selected(pair1, matched_intensity_references)

    pair_address_has_invalid_intensity_reference = (self.params.scaling.i_model.sigmas() < 0.).select(pair0)
    pair_addresses_with_invalid_intensity_reference = pair_address_has_invalid_intensity_reference.iselection(True)
    reflection_addresses_with_invalid_intensity_reference = pair1.select(pair_addresses_with_invalid_intensity_reference)

    use_weights = False # New facility for getting variance-weighted correlation

    if use_weights:
      # variance weighting
      # obviously there is a chance of divide by zero here so recheck the code
      intensity_weight = 1./reflections['intensity.sum.variance'].select(pair1)
    else:
      IW = intensity_weight = flex.double(len(reflections), 1.)

    intensity_weight.set_selected(reflection_addresses_with_invalid_intensity_reference, 0.)
    intensity_weight.set_selected(new_matches.singles(1), 0.) # observations with no corresponding reference intensity

    reflections["intensity_weight"] = intensity_weight
    reflections["intensity_reference"] = intensity_reference
    reflections["resolution_d"] = miller_set.unit_cell().d(reflections['miller_index_asymmetric']) # Have to choose a unique unit cell; despite unit cell variation.
    reflections["resolution_DSSQ"] = miller_set.unit_cell().d_star_sq(reflections['miller_index_asymmetric'])

    # Now get the base crystal orientations and beam properties
    reflections["expt_idx"] = flex.size_t([experiments.find(reflections[refl_id]['exp_id']) for refl_id in range(len(reflections))])
    reflections["A_matrix"] = flex.mat3_double( [C.get_A() for C in experiments.crystals()] ).select(reflections["expt_idx"])
    reflections["s0_vec"] = flex.vec3_double( [B.get_s0() for B in experiments.beams()] ).select(reflections["expt_idx"])
    reflections["wavelength"] = flex.double( [B.get_wavelength() for B in experiments.beams()] ).select(reflections["expt_idx"])


    # end of reference column deep dive, resuling in new columns for the reflection table

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI
    incoming_experiment_count = comm.reduce(len(experiments), MPI.SUM, 0)
    incoming_counts = comm.gather(len(experiments), root=0)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Incoming to postrefine with %d experiments"%(incoming_experiment_count))
      incoming_counts = flex.double(incoming_counts)
      self.logger.main_log("Min: %d  Max: %d  Mean: %.1f"%(flex.min(incoming_counts),flex.max (incoming_counts),flex.mean(incoming_counts)))

    incoming_reflection_count = comm.reduce(len(reflections), MPI.SUM, 0)
    incoming_refls = comm.gather(len(reflections), root=0)
    if self.mpi_helper.rank == 0:
      self.logger.main_log("Incoming to postrefine with %d reflections"%(incoming_reflection_count))
      incoming_refls = flex.double(incoming_refls)
      self.logger.main_log("Min: %d  Max: %d  Mean: %.1f"%(flex.min(incoming_refls),flex.max (incoming_refls),flex.mean(incoming_refls)))

    # parameters for each experiment: G, BFACTOR.  2 parameters x # of experiments
    # common parameters over all experiments: eta, deff.  + 2 common
    mpi_refinery = ratchet1_refinery.from_specific_and_common(
                     specific=["G","BFACTOR"], common=["eta","deff"],
                     experiments=experiments, reflections=reflections,
                     mpi_helper = self.mpi_helper,
                   )
    mpi_refinery.logger = self.logger
    mpi_refinery.compute_starting_values(self.mpi_helper,self.params)
    mpi_refinery.run_ratchet(self.mpi_helper,self.params)
    # BEGIN LEGACY for loop through experiments. start here.  create a master root=0 minimizer

    for experiment in experiments:

      exp_reflections = reflections.select(reflections['exp_id'] == experiment.identifier)
      notional_G_value = exp_reflections['G_factor'][0]
      notional_B_value = exp_reflections['BFACTOR'][0]

      # {'miller_index_asymmetric': (28, 12, 13), 'exp_id': 'e9d5627f1140e5eb92c036e3ef6ca78b', 'intensity.sum.value': -13.592601855291326, 'miller_index': (-28, 12, 13), 'intensity.sum.variance': 1690.3049948003838}

      # Build a miller array for the experiment reflections with original miller indexes
      exp_miller_indices_original = miller.set(target_symm, exp_reflections['miller_index'], not self.params.merging.merge_anomalous)
      observations_original_index = miller.array(exp_miller_indices_original,
                                                 exp_reflections['intensity.sum.value'],
                                                 flex.double(flex.sqrt(exp_reflections['intensity.sum.variance'])))

      assert exp_reflections.size() == exp_miller_indices_original.size()
      assert observations_original_index.size() == exp_miller_indices_original.size()

      # Build a miller array for the experiment reflections with asu miller indexes
      exp_miller_indices_asu = miller.set(target_symm, exp_reflections['miller_index_asymmetric'], True)
      observations = miller.array(exp_miller_indices_asu, exp_reflections['intensity.sum.value'], flex.double(flex.sqrt(exp_reflections['intensity.sum.variance'])))

      matches = miller.match_multi_indices(miller_indices_unique = miller_set.indices(), miller_indices = observations.indices())

      pair1 = flex.int([pair[1] for pair in matches.pairs()]) # refers to the observations
      pair0 = flex.int([pair[0] for pair in matches.pairs()]) # refers to the model

      assert exp_reflections.size() == exp_miller_indices_original.size()
      assert observations_original_index.size() == exp_miller_indices_original.size()

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
      """

      # RB: By design, for MPI-Merge "include negatives" is implicitly True
      SWC = simple_weighted_correlation(I_weight, I_reference, I_observed)
      if self.params.output.log_level == 1:
        self.logger.log("Old correlation is: %f"%SWC.corr)

      if self.params.postrefinement.algorithm == "ratchet1":

        Rhall = flex.double()

        for mill in MILLER:
          H = matrix.col(mill)
          Xhkl = Astar*H
          Rh = ( Xhkl + BEAM ).length() - (1./WAVE)
          Rhall.append(Rh)

        Rs = math.sqrt(flex.mean(Rhall*Rhall))

        RS = 1./400. # reciprocal effective domain size of 1 micron
        #RS = Rs        # try this empirically determined approximate, monochrome, a-mosaic value
        if self.mpi_helper.rank in range(20):
          for MAN in mpi_refinery.parameter_managers:
            if MAN.exp_id != experiment.identifier: continue
            if MAN.type != "G": continue
            print([SWC.slope, BFACTOR, RS, 0., 0.], "of", len(experiments),MAN.value)
            assert approx_equal(SWC.slope,MAN.initial_value_actual)



        current = flex.double([SWC.slope, BFACTOR, RS, 0., 0.])

        parameterization_class = rs_parameterization
        refinery = rs_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE, ICALCVEC = I_reference, IOBSVEC = I_observed)

      func = refinery.fvec_callable(parameterization_class(current))
      functional = flex.sum(func * func)

      if self.params.output.log_level == 1:
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
        print ("Experiment %s, Global G=%8.5f Local G=%8.5f, Global B=%8.5f Local B=%8.5f"%(
               experiment.identifier, notional_G_value, self.MINI.x[0], notional_B_value, self.MINI.x[1]
              ))

      except (AssertionError, ValueError, RuntimeError) as e:
        error_detected = True
        reason = repr(e)
        if not reason:
          reason = "Unknown error"
        if not reason in experiments_rejected_by_reason:
          experiments_rejected_by_reason[reason] = 1
        else:
          experiments_rejected_by_reason[reason] += 1

      if not error_detected:
        new_experiments.append(experiment)

        new_exp_reflections = flex.reflection_table()
        new_exp_reflections['miller_index_asymmetric']  = flex.miller_index(result_observations.indices())
        new_exp_reflections['intensity.sum.value']      = flex.double(result_observations.data())
        new_exp_reflections['intensity.sum.variance']   = flex.double(flex.pow(result_observations.sigmas(),2))
        new_exp_reflections['exp_id']                   = flex.std_string(len(new_exp_reflections), experiment.identifier)
        new_reflections.extend(new_exp_reflections)
      '''
      # debugging
      elif reason.startswith("ValueError"):
        self.logger.log("Rejected b/c of value error exp id: %s; unit cell: %s"%(exp_id, str(experiment.crystal.get_unit_cell())) )
      '''

    # END of LEGACY for loop, through experiments

    # report rejected experiments, reflections
    experiments_rejected_by_postrefinement = len(experiments) - len(new_experiments)
    reflections_rejected_by_postrefinement = reflections.size() - new_reflections.size()

    self.logger.log("Experiments rejected by post-refinement: %d"%experiments_rejected_by_postrefinement)
    self.logger.log("Reflections rejected by post-refinement: %d"%reflections_rejected_by_postrefinement)

    all_reasons = []
    for reason, count in six.iteritems(experiments_rejected_by_reason):
      self.logger.log("Experiments rejected due to %s: %d"%(reason,count))
      all_reasons.append(reason)

    comm = self.mpi_helper.comm
    MPI = self.mpi_helper.MPI

    # Collect all rejection reasons from all ranks. Use allreduce to let each rank have all reasons.
    all_reasons  = comm.allreduce(all_reasons, MPI.SUM)
    all_reasons = set(all_reasons)

    # Now that each rank has all reasons from all ranks, we can treat the reasons in a uniform way.
    total_experiments_rejected_by_reason = {}
    for reason in all_reasons:
      rejected_experiment_count = 0
      if reason in experiments_rejected_by_reason:
        rejected_experiment_count = experiments_rejected_by_reason[reason]
      total_experiments_rejected_by_reason[reason] = comm.reduce(rejected_experiment_count, MPI.SUM, 0)

    total_accepted_experiment_count = comm.reduce(len(new_experiments), MPI.SUM, 0)

    # how many reflections have we rejected due to post-refinement?
    rejected_reflections = len(reflections) - len(new_reflections);
    total_rejected_reflections = self.mpi_helper.sum(rejected_reflections)

    if self.mpi_helper.rank == 0:
      for reason, count in six.iteritems(total_experiments_rejected_by_reason):
        self.logger.main_log("Total experiments rejected due to %s: %d"%(reason,count))
      self.logger.main_log("Total experiments accepted: %d"%total_accepted_experiment_count)
      self.logger.main_log("Total reflections rejected due to post-refinement: %d"%total_rejected_reflections)

    self.logger.log_step_time("POSTREFINEMENT", True)

    return new_experiments, new_reflections

  def run_plain(self):

    out = StringIO()
    self.MINI = lbfgs_minimizer_base(self.params,
                                     current_x = self.current,
                                     parameterization = self.parameterization_class,
                                     refinery = self.refinery,
                                     out = out)
    if self.params.output.log_level == 1:
      self.logger.log("\n" + out.getvalue())

  def result_for_cxi_merge(self):

    scaler = self.refinery.scaler_callable(self.parameterization_class(self.MINI.x))

    if self.params.postrefinement.algorithm == "ratchet1":
      fat_selection = (self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x)) > 0.2)
      fats = self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x))
    else:
      fat_selection = (self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x)) < 0.9)

    fat_count = fat_selection.count(True)

    # reject an experiment with insufficient number of near-full reflections
    if fat_count < 3:

      if self.params.output.log_level == 1:
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

    if self.params.output.log_level == 1:
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

class parameter_manager(group_args): # abstract class
    def __init__(self, **kwargs):
      group_args.__init__(self,**kwargs)
      mandatory = ["exp_id", "param_id"]
      for key in mandatory: getattr(self,key)
      self.refl_selection_already_set = False

    def set_from_current(self, current):
      self.value = current[self.param_id]
    def set_refl_selection(self, refls):
      if self.scope == "common":
        self.refl_selection_already_set = flex.bool(len(refls["exp_id"]),True)
        return
      self.refl_selection_already_set = (refls["exp_id"]==self.exp_id)
    def get_jacobian_column(self):
      value = self.get_jacobian_column_refined()
      if value == "Not implemented":
        return flex.double(len(self.reflections)) # default zero, no gradient.
      else:
        return value
    def get_jacobian_column_refined(self):
      return "Not implemented"

class G_manager(parameter_manager):
    def __init__(self,**kwargs):
      parameter_manager.__init__(self,**kwargs)
      self.type = "G"
      self.scope = "specific"
    def initial_value(self,params,mpi_helper,**kwargs):
      exp_reflections = self.reflections.select(self.reflections['exp_id'] == self.exp_id)
      alternate_SWC = simple_weighted_correlation(exp_reflections["intensity_weight"],
                                                  exp_reflections["intensity_reference"],
                                                  exp_reflections["intensity.sum.value"])
      self.value = alternate_SWC.slope
      self.initial_value_actual = alternate_SWC.slope # for debug only
      return self.value
    def get_jacobian_column_refined(self):
      Z = flex.double(len(self.reflections)).set_selected(self.refl_selection_already_set, self.reflections["G_terms"])
      return Z

class BFACTOR_manager(parameter_manager):
    def __init__(self,**kwargs):
      parameter_manager.__init__(self,**kwargs)
      self.type = "BFACTOR"
      self.scope = "specific"
    def initial_value(self,**kwargs):
      self.value = 0.
      return self.value # resolution-dependent intensity falloff in Angstrom**2
    def get_jacobian_column_refined(self):
      Z = flex.double(len(self.reflections)).set_selected(self.refl_selection_already_set, self.reflections["B_terms"])
      return Z

class deff_manager(parameter_manager):
    def __init__(self,**kwargs):
      parameter_manager.__init__(self,**kwargs)
      self.type = "deff"
      self.scope = "common"
    def initial_value(self,**kwargs):
      self.value = 400.
      return self.value # effective domain size in Angstroms
    def get_jacobian_column_refined(self):
      Z = self.reflections["P_terms"] * self.reflections["d_PART_d_rs"] * (-1. / (self.value * self.value) )
      return Z

class eta_manager(parameter_manager):
    def __init__(self,**kwargs):
      parameter_manager.__init__(self,**kwargs)
      self.type = "eta"
      self.scope = "common"
    def initial_value(self,**kwargs):
      self.value = 0.
      return self.value # effective mosaic angle in radians
    def get_jacobian_column_refined(self):
      Z = self.reflections["P_terms"] * self.reflections["d_PART_d_rs"] * ( 0.5 / self.reflections["resolution_d"] )
      return Z

manager_type = {"G":G_manager, "BFACTOR":BFACTOR_manager, "eta":eta_manager, "deff":deff_manager}

class ratchet1_refinery(refinery_base):

    def __init__(self, **kwargs):
      group_args.__init__(self,**kwargs)
      # ratchet1_refinery operates over many-experiments, so requires slightly different
      #   mandatory inputs than single-experiment refinery.
      #   New single single columns will be required; still yet to be worked out
      mandatory = [] # "ORI","MILLER","BEAM","WAVE","ICALCVEC","IOBSVEC"] # example only, will figure out
      for key in mandatory: getattr(self,key)
      self.lookup = {}
      for idx,item in enumerate(self.parameter_managers):
        self.lookup[getattr(item,"exp_id","")+item.type] = idx

    def compute_starting_values(self, mpi_helper, params):
      comm = mpi_helper.comm
      this_rank_initial_values = flex.double()
      for manager in self.parameter_managers:
        if manager.scope=="specific":  this_rank_initial_values.append(manager.initial_value(mpi_helper=mpi_helper,params=params))
      reports = comm.gather(this_rank_initial_values, root = 0)
      current = flex.double()
      if mpi_helper.rank==0:
        for rank_report in reports:  current = current.concatenate(rank_report) # better way?
        for manager in self.parameter_managers:
          if manager.scope=="common":  current.append(manager.initial_value())
      self.starting_x = comm.bcast(current, root=0)
      #sets the initial starting values

    def run_ratchet(self,mpi_helper, params):

      out = StringIO()
      self.MINI = lbfgs_minimizer_ratchet(params,
                                          current_x = self.starting_x,
                                          parameterization = None,
                                          refinery_target = self,
                                          out = out)
      if params.output.log_level == 1:
        self.logger.log("\n" + out.getvalue())


    @classmethod
    def from_specific_and_common(cls, specific, common, experiments, reflections, mpi_helper):
      comm = mpi_helper.comm
      param_counts = flex.int(comm.allgather(len(specific)*len(experiments)))
      total_specific = flex.sum(param_counts)
      rank_start = flex.sum(param_counts[:mpi_helper.rank])
      print ("Rank %d starts at %d of %d"%(mpi_helper.rank, rank_start, total_specific))
      parameter_managers = []
      for experiment in experiments:
        for param_type in specific:
          access_hash = "%s_%s"%(experiment.identifier, param_type) # not used, delete later
          parameter_managers.append( manager_type[param_type](
            exp_id = experiment.identifier,
            param_id = rank_start,
            experiments = experiments,
            reflections = reflections
          ))
          rank_start += 1
      for param_type in common:
          parameter_managers.append( manager_type[param_type](
            exp_id = "",
            param_id = total_specific,
            experiments = experiments,
            reflections = reflections
          ))
          total_specific+=1
      # By the way, it is a bad idea to use MD5 hashes as experiment identifiers
      # due to MD5 collision vulnerability.  What is the probability that 1000 users
      # each processing a 1-million image data set will all together experience
      # 1 collision? Why not just use timestamp that never overlaps?
      return cls(parameter_managers=parameter_managers,
                 experiments=experiments, reflections=reflections,
                 total_specific = len(specific)*len(experiments))

    def take_current_values(self,current):
      for managed_parameter in self.parameter_managers:
        managed_parameter.set_from_current(current)

    def fvec_callable(self):
      # in MPI design, values are already transferred to the parameter_managers in take_current_values()
      # next step is to populate the flex vectors needed to calculate target fvec.
      # This needs to be done every time the parameter set is updated (take_current_values() invalidates these vectors):
      self.reflections["G_factor"] = flex.double([self.parameter_managers[self.lookup[e+"G"]].value for e in self.reflections["exp_id"]])
      self.reflections["BFACTOR"] = flex.double([self.parameter_managers[self.lookup[e+"BFACTOR"]].value for e in self.reflections["exp_id"]])

      self.RLP_model = "PG" # PB="Lorentzian"  PG="Gaussian"
      PART = self.get_partiality_array(RLP_model = self.RLP_model, refls=self.reflections)
      EXP = flex.exp(-2.* self.reflections["BFACTOR"] * self.reflections["resolution_DSSQ"])
      residuals = (self.reflections["G_factor"] * EXP * PART * self.reflections["intensity_reference"] - self.reflections["intensity.sum.value"])

      return residuals

    def get_partiality_array(self,RLP_model, refls):
      # rs = 1/Deff + eta/(2d)
      # not sure how to design parameter access
      Deff = [p for p in self.parameter_managers if p.type=="deff"][0].value
      eta = [p for p in self.parameter_managers if p.type=="eta"][0].value
      # vector formula, Sauter (2015) equation (6):
      rs = (1./refls["resolution_d"]) * (eta/2.) + (1./Deff)
      rh = self.get_Rh_array(refls)
      rs_sq = rs*rs
      rh_sq = rh*rh

      # best design choice to make these arrays permanent in the reflection table?
      # or should they become class attributes?
      self.rs = rs
      self.rh = rh
      self.rs_sq = rs_sq
      self.rh_sq = rh_sq

      if RLP_model=="PB":
        return rs_sq / ( (2. * rh_sq) + rs_sq )
      else:
        assert RLP_model=="PG"
        return flex.exp( -2.*math.log(2.) * rh_sq / rs_sq )

    def get_Rh_array(self,refls):
      eff_Astar = self.get_eff_Astar(refls)
      h = refls["miller_index"].as_vec3_double()
      x = eff_Astar * h
      Svec = x + refls["s0_vec"]
      Rh = Svec.norms() - (1./refls["wavelength"])
      return Rh

    def get_eff_Astar(self,refls):
      # thetax and thetay are no-ops since orientation is fixed
      return refls["A_matrix"]

    def jacobian_callable(self):
      for managed_parameter in self.parameter_managers:
        if not managed_parameter.refl_selection_already_set:
          managed_parameter.set_refl_selection(self.reflections) # sets which reflections determine a given parameter

      PART = self.get_partiality_array(RLP_model=self.RLP_model, refls=self.reflections)
      EXP = flex.exp(-2. * self.reflections["BFACTOR"] * self.reflections["resolution_DSSQ"])
      self.reflections["G_terms"] = (EXP * PART * self.reflections["intensity_reference"])
      self.reflections["B_terms"] = (self.reflections["G_factor"] * EXP * PART * self.reflections["intensity_reference"])*(
                                     -2.* self.reflections["resolution_DSSQ"] )
      self.reflections["P_terms"] = (self.reflections["G_factor"] * EXP * self.reflections["intensity_reference"])

      if self.RLP_model=="PG":
        self.reflections["d_PART_d_rs"] = PART * (4. * math.log(2.)) * self.rh_sq / (self.rs_sq * self.rs)
      else:
        denom_factor = 2. * self.rh_sq + self.rs_sq
        self.reflections["d_PART_d_rs"] = 4. * self.rh_sq * self.rs / ( denom_factor * denom_factor )

      jacobian_columns = []
      for managed_parameter in self.parameter_managers:
        jacobian_columns.append(managed_parameter.get_jacobian_column())

      #P_terms = (values.G * EXP * self.ICALCVEC)

      #thetax = values.thetax; thetay = values.thetay;
      #Rx = matrix.col((1,0,0)).axis_and_angle_as_r3_rotation_matrix(thetax)
      #dRx_dthetax = matrix.col((1,0,0)).axis_and_angle_as_r3_derivative_wrt_angle(thetax)
      #Ry = matrix.col((0,1,0)).axis_and_angle_as_r3_rotation_matrix(thetay)
      #dRy_dthetay = matrix.col((0,1,0)).axis_and_angle_as_r3_derivative_wrt_angle(thetay)
      #ref_ori = matrix.sqr(self.ORI.reciprocal_matrix())
      #miller_vec = self.MILLER.as_vec3_double()
      #ds1_dthetax = flex.mat3_double(len(self.MILLER),Ry * dRx_dthetax * ref_ori) * miller_vec
      #ds1_dthetay = flex.mat3_double(len(self.MILLER),dRy_dthetay * Rx * ref_ori) * miller_vec

      #s1vec = self.get_s1_array(values)
      #s1lenvec = flex.sqrt(s1vec.dot(s1vec))
      #dRh_dthetax = s1vec.dot(ds1_dthetax)/s1lenvec
      #dRh_dthetay = s1vec.dot(ds1_dthetay)/s1lenvec
      #rs = values.RS
      #Rh = self.get_Rh_array(values)
      #rs_sq = rs*rs
      #denomin = (2. * Rh * Rh + rs_sq)
      #dPB_dRh = { "lorentzian": -PB * 4. * Rh / denomin,
      #            "gaussian": -PB * 4. * math.log(2) * Rh / rs_sq }[self.profile_shape]
      #dPB_dthetax = dPB_dRh * dRh_dthetax
      #dPB_dthetay = dPB_dRh * dRh_dthetay
      #Px_terms = P_terms * dPB_dthetax; Py_terms = P_terms * dPB_dthetay

      return jacobian_columns # [ G, B, G, B, ..., eta, deff ]


from xfel.merging.application.postrefine.postrefinement import unpack_base

class rs_parameterization(unpack_base):
  def __init__(YY,values):
    super(rs_parameterization, YY).__init__(values)
    YY.keys = ['G', 'BFACTOR','RS','thetax','thetay']

  def show(YY, out):
    print ("G: %10.7f; B: %10.7f; RS: %10.7f; THETAX: %7.3f deg; THETAY: %7.3f deg"\
          %(YY.G, YY.BFACTOR, YY.RS, 180.*YY.thetax/math.pi, 180.*YY.thetay/math.pi)\
          , file=out)

class lbfgs_minimizer_ratchet(object):

  def __init__(self, params, current_x=None, parameterization=None, refinery_target=None, out=None,
               min_iterations=0, max_calls=50, max_drop_eps=1.e-5):  # reset to 100 max calls after debug
    adopt_init_args(self, locals())
    self.n = current_x.size()
    self.x = current_x
    from scitbx import lbfgs
    from scitbx.lbfgs.tst_mpi_split_evaluator import mpi_split_evaluator_run
    P = Profiler("LBFGS total run")
    self.minimizer = mpi_split_evaluator_run(
      target_evaluator=self,
      termination_params=lbfgs.termination_parameters(
        traditional_convergence_test=True,
        traditional_convergence_test_eps=1.e-1,
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
    P = Profiler("functional and gradients only")
    # take current parameters self.x, translate into the data structure of the refinery_target
    self.refinery_target.take_current_values(self.x)
    self.func = self.refinery_target.fvec_callable()
    functional = flex.sum(self.refinery_target.reflections[ "intensity_weight" ] *self.func*self.func)
    self.f = functional
    jacobian = self.refinery_target.jacobian_callable()
    self.g = flex.double(self.n)
    for ix in range(len(self.refinery_target.parameter_managers)):
      self.g[ self.refinery_target.parameter_managers[ix].param_id ] = flex.sum(
         2. * self.refinery_target.reflections[ "intensity_weight" ] * self.func * jacobian[ix])
    print("rms %10.3f"%math.sqrt(flex.sum(self.refinery_target.reflections[ "intensity_weight" ]*self.func*self.func)/
                                 flex.sum(self.refinery_target.reflections[ "intensity_weight" ])), end=' ')#, file=self.out)
    print(list(self.x[-4:]))#,file=self.out)
    return self.f, self.g


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

    if self.params.postrefinement.algorithm == 'ratchet1':
      for p in ['RS','thetax','thetay']:# self.params.postrefinement.rs.fix:
      #for p in ['RS','BFACTOR','thetax','thetay']:# self.params.postrefinement.rs.fix:
        #print ('fixed-->',p,values.keys,values.keys.index(p))
        self.g[values.keys.index(p)] = 0

    print ("rms %10.3f; "%math.sqrt(flex.mean(self.func*self.func)), file=self.out, end='')
    values.show(self.out)

    return self.f, self.g
  """start here
  DONE put some instrumentation in the postrefinement
  DONE gather all parameters, and initialize them all
  NOT MODULAR, HOLD OFF FOR NOW clean out the ratchet file so it imports most stuff
  create a master root=0 minimizer
    CHECKED first check that I_model columns can be created with proper invalid flags
    OK FOR NOW check that minimizer can work refining only G,B
  sept 22 path
    DONE implement default of PG instead of PB
    DONE double check analytical derivatives for G, B
    DONE test if refined G values match the old refinery.  (under PB model)
    DONE same for B
      result:  while the 1-image vs. 200-image tests did not give equal parameters (G==G, B==B)
      the results were highly correlated, > 95%.  Seems good enough.
    DONE WITH 37 RANKS scale up to full dataset & test performance
      result:  50 seconds in LBFGS, 1125 seconds f,g in longest rank, 20 seconds in shortest
      10x improvement if load is rebalanced
      10x to 50x if recoded in C++
    DONE write derivatives for ETA, DEFF
      Eta, Deff do not refine to realistic value
      hypotheses:
        derivatives are wrong
        math is wrong somewhere else
        eta derivative should be logarithmic so it remains positive
        Deff derivative should be log so it scales correctly
        Use LSQ curvature approximation so all parameters scale correctly.  Can it be made to work with MPI?
  new features desired:
    refinement of DEFF
    refinement of ETA
    weighted treatment of observations
    choice of weighting
    choice of Gauss / Lorentz
    analytical derivatives
  new rejection regime:
    some images rejected because minimization problem can't be set up
    some rejected because criteria fail after minimization (fat count, eg)
    some raise exception during refinement (try to engineer this category away)

  """
  """
Global postrefinement prototype

Experiment-specific G,B along with global single values of eta, Deff
Per-experiment parameters refined all together in single LBFGS minimizer, no sequential loop
Use MPI to gather per-experiment contributions.
Known bug: cannot use more MPI ranks than input expt/refl files
Use case of 5874 images; takes 2000 sec for scitbx.flex-based vector ops for funtion & analytical gradient
Only 50 seconds for LBFGS engine (100 steps, not converged)
Eta, Deff do not refine; curvatures not yet implemented but might be needed
Variance weighting stub is present, but needs to be implemented next
Choice between Lorentzian or Gaussian RLP model
Ability to disable any parameter refinement by setting column return value to default
Refined parameters not yet applied and results not yet output, stills needs to be done
  """

  def __del__(self):
    values = self.parameterization(self.x)
    print ("FINALMODEL", file=self.out)
    print ("rms %10.3f; "%math.sqrt(flex.mean(self.func*self.func)), file=self.out, end='')
    values.show(self.out)

if __name__ == '__main__':
  from xfel.merging.application.worker import exercise_worker
  exercise_worker(ratchet1_postrefinement)
