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
from six.moves import cStringIO as StringIO
from libtbx.development.timers import Profiler

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

    use_weights = True # New facility for getting variance-weighted correlation

    if use_weights:
      # Weighting schemes
      # obviously there is a chance of divide by zero here so recheck the code
      intensity_weight = 1./reflections['intensity.sum.variance']#.select(pair1) # Variance
      #intensity_weight = flex.abs(reflections['intensity.sum.value'])/reflections['intensity.sum.variance'] # Gentle
      #intensity_weight = flex.pow(reflections['intensity.sum.value'],2)/reflections['intensity.sum.variance'] # Extreme
    else:
      intensity_weight = flex.double(len(reflections), 1.)

    intensity_weight.set_selected(reflection_addresses_with_invalid_intensity_reference, 0.)
    intensity_weight.set_selected(new_matches.singles(1), 0.) # observations with no corresponding reference intensity

    reflections["intensity_weight"] = intensity_weight
    reflections["intensity_reference"] = intensity_reference
    reflections["resolution_d"] = miller_set.unit_cell().d(reflections['miller_index_asymmetric']) # Have to choose a unique unit cell; despite unit cell variation.
    reflections["resolution_DSSQ"] = miller_set.unit_cell().d_star_sq(reflections['miller_index_asymmetric'])

    # Now get the base crystal orientations and beam properties
    reflections["expt_idx"] = flex.size_t([experiments.find(reflections[refl_id]['exp_id']) for refl_id in range(len(reflections))])
    reflections["A_matrix"] = flex.mat3_double( [C.get_A() for C in experiments.crystals()] ).select(reflections["expt_idx"])
    reflections["s0_vec"] = flex.vec3_double( [e.beam.get_s0() for e in experiments] ).select(reflections["expt_idx"])
    reflections["wavelength"] = flex.double( [e.beam.get_wavelength() for e in experiments] ).select(reflections["expt_idx"])

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
    #mpi_refinery.run_ratchet(self.mpi_helper,self.params)
    mpi_refinery.run_ratchet_levmar(self.mpi_helper,self.params)

    new_experiments, new_reflections = self.apply_new_parameters(mpi_refinery)

    # report rejected experiments, reflections
    experiments_rejected_by_postrefinement = len(experiments) - len(new_experiments)
    reflections_rejected_by_postrefinement = reflections.size() - new_reflections.size()

    self.logger.log("Experiments rejected by post-refinement: %d"%experiments_rejected_by_postrefinement)
    self.logger.log("Reflections rejected by post-refinement: %d"%reflections_rejected_by_postrefinement)

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
    if self.params.output.log_level == 1:
      self.logger.log("\n" + out.getvalue())

  def apply_new_parameters(self, refinery_target):
    reflections = refinery_target.reflections

    corrections = refinery_target.get_corrections()
    bad_corrections = corrections <= 1e-2
    self.logger.log("Removing %d reflections with unstable intensity corrections out of %d"%(bad_corrections.count(True), len(reflections)))
    refinery_target.reflections = reflections = reflections.select(~bad_corrections)
    corrections = corrections.select(~bad_corrections)

    partiality = refinery_target.get_partiality_array(RLP_model=refinery_target.RLP_model,
                                                      refls=reflections)

    reflections['intensity.sum.value'] /= corrections
    sigmas = flex.sqrt(reflections['intensity.sum.variance'])
    reflections['intensity.sum.variance'] = (sigmas/corrections)**2

    assert self.params.postrefinement.algorithm == "ratchet1"
    fat_selection = partiality > 0.2

    new_experiments = ExperimentList()
    new_reflections = flex.reflection_table()
    param_deff = refinery_target.parameter_managers[-1]
    param_eta = refinery_target.parameter_managers[-2]

    for expt_idx, experiment in enumerate(refinery_target.experiments):
      sel = reflections['expt_idx'] == expt_idx
      n_refl = sel.count(True)
      fat_count = (fat_selection & sel).count(True)

      # reject an experiment with insufficient number of near-full reflections
      if fat_count < 3:
        if self.params.output.log_level == 1:
          self.logger.log("Rejected experiment %d, because: On total %5d the fat selection is %5d"%(expt_idx, n_refl, fat_count))
          continue

      if self.params.output.log_level == 1:
        self.logger.log("Expt %d: On total %5d the fat selection is %5d"%(expt_idx, n_refl, fat_count))

      new_experiments.append(experiment)
      refls = reflections.select(fat_selection & sel)
      refls['expt_idx'] = flex.size_t(len(refls), len(new_experiments)-1)
      new_reflections.extend(refls)

      experiment.crystal.set_domain_size_ang(2/math.exp(param_deff.value))
      experiment.crystal.set_half_mosaicity_deg(180/math.pi * 0.5 * math.exp(param_eta.value))
    return new_experiments, new_reflections

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
      self.value = 100.
      return self.value # effective domain size in Angstroms
    def get_jacobian_column_refined(self):
      Z = self.reflections["P_terms"] * self.reflections["d_PART_d_rs"] * (-1. / (self.value * self.value) )
      return Z

class log_alpha_manager(deff_manager):
    def initial_value(self,**kwargs):
      super(log_alpha_manager, self).initial_value(**kwargs)
      self.value = math.log(2/self.value) # u = ln(alpha) = ln(2/Deff)
      return self.value # effective domain size in Angstroms
    def get_jacobian_column_refined(self):
      Z = self.reflections["P_terms"] * self.reflections["d_PART_d_rs"] * (0.5*math.exp(self.value))
      return Z

class eta_manager(parameter_manager):
    def __init__(self,**kwargs):
      parameter_manager.__init__(self,**kwargs)
      self.type = "eta"
      self.scope = "common"
    def initial_value(self,**kwargs):
      self.value = 0.2
      return self.value # effective mosaic angle in radians
    def get_jacobian_column_refined(self):
      Z = self.reflections["P_terms"] * self.reflections["d_PART_d_rs"] * ( 0.5 / self.reflections["resolution_d"] )
      return Z

class log_eta_manager(eta_manager):
    def initial_value(self,**kwargs):
      super(log_eta_manager, self).initial_value(**kwargs)
      self.value = math.log(self.value) # u = ln(eta)
      return self.value # effective mosaic angle in radians
    def get_jacobian_column_refined(self):
      Z = self.reflections["P_terms"] * self.reflections["d_PART_d_rs"] * ( math.exp(self.value) * 0.5 / self.reflections["resolution_d"] )
      return Z

manager_type = {"G":G_manager, "BFACTOR":BFACTOR_manager, "eta":log_eta_manager, "deff":log_alpha_manager}

class ratchet1_refinery(group_args):
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

    def run_ratchet_levmar(self,mpi_helper, params):
      out = StringIO()
      self.MINI = levmar_minimizer_ratchet(current_x = self.starting_x,
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

    def get_corrections(self):
      self.RLP_model = "PG" # PB="Lorentzian"  PG="Gaussian"
      PART = self.get_partiality_array(RLP_model = self.RLP_model, refls=self.reflections)
      EXP = flex.exp(-2.* self.reflections["BFACTOR"] * self.reflections["resolution_DSSQ"])
      return self.reflections["G_factor"] * EXP * PART

    def fvec_callable(self):
      # in MPI design, values are already transferred to the parameter_managers in take_current_values()
      # next step is to populate the flex vectors needed to calculate target fvec.
      # This needs to be done every time the parameter set is updated (take_current_values() invalidates these vectors):
      self.reflections["G_factor"] = flex.double([self.parameter_managers[self.lookup[e+"G"]].value for e in self.reflections["exp_id"]])
      self.reflections["BFACTOR"] = flex.double([self.parameter_managers[self.lookup[e+"BFACTOR"]].value for e in self.reflections["exp_id"]])

      residuals = (self.get_corrections() * self.reflections["intensity_reference"] - self.reflections["intensity.sum.value"])

      return residuals

    def get_partiality_array(self,RLP_model, refls):
      # rs = 1/Deff + eta/(2d)
      # not sure how to design parameter access
      u = [p for p in self.parameter_managers if p.type=="deff"][0].value
      #Deff = [p for p in self.parameter_managers if p.type=="deff"][0].value
      alpha = math.exp(u)
      try:
        u = [p for p in self.parameter_managers if p.type=="eta"][0].value
        #eta = [p for p in self.parameter_managers if p.type=="eta"][0].value
      except IndexError:
        u = 0
        #eta = 0
      eta = math.exp(u)

      # vector formula, Sauter (2015) equation (6):
      rs = (1./refls["resolution_d"]) * (eta/2.) + (alpha/2)
      #rs = (1./refls["resolution_d"]) * (eta/2.) + (1./Deff)
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
    #self.diag_mode = "always"
    self.outliers = False
    from scitbx import lbfgs
    from scitbx.lbfgs.tst_mpi_split_evaluator import mpi_split_evaluator_run
    self.diag_mode=None # None, "once", or "always"
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

  def do_outliers(self):
    if not self.outliers:
      from scitbx.math import five_number_summary
      sel = self.refinery_target.reflections[ "intensity_weight" ] != 0
      print ("Initial count rejected", sel.count(False), "out of", len(sel))
      min_x, q1_x, med_x, q3_x, max_x = five_number_summary(self.func.select(sel))
      print ("Five number summary: min %.1f, q1 %.1f, med %.1f, q3 %.1f, max %.1f"%(min_x, q1_x, med_x, q3_x, max_x))
      iqr_x = q3_x - q1_x
      cut_x = 30 * iqr_x
      sel.set_selected(self.func > q3_x + cut_x, False)
      sel.set_selected(self.func < q1_x - cut_x, False)
      self.refinery_target.reflections[ "intensity_weight" ].set_selected(~sel, 0)
      print ("Rejecting", sel.count(False), "out of", len(sel))
      self.outliers = True

  def compute_functional_and_gradients(self):
    P = Profiler("functional and gradients only")
    # take current parameters self.x, translate into the data structure of the refinery_target
    self.refinery_target.take_current_values(self.x)
    self.func = self.refinery_target.fvec_callable()
    self.do_outliers()

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


    # check finite differences
    """
    DELTA = 1.e-9
    ran = list(range(len(self.refinery_target.parameter_managers)))
    for param_id in ran[0:6] + ran[-10:]:
      parameter = self.refinery_target.parameter_managers[param_id]
      current = parameter.value
      parameter.value += DELTA
      func = self.refinery_target.fvec_callable()
      dfunctional = flex.sum(self.refinery_target.reflections[ "intensity_weight" ]*func*func)
      print ("%10s"%parameter.type, param_id, functional, dfunctional, self.g[param_id], (dfunctional-functional)/DELTA)
      parameter.value = current
    """

    return self.f, self.g

  def compute_functional_gradients_diag(self):
    P = Profiler("functional, gradients and diagonal")
    # take current parameters self.x, translate into the data structure of the refinery_target
    self.refinery_target.take_current_values(self.x)
    self.func = self.refinery_target.fvec_callable()
    self.do_outliers()
    functional = flex.sum(self.refinery_target.reflections[ "intensity_weight" ] *self.func*self.func)
    self.f = functional
    jacobian = self.refinery_target.jacobian_callable()
    self.g = flex.double(self.n)
    self.c = flex.double(self.n)
    for ix in range(len(self.refinery_target.parameter_managers)):
      d = 2. * self.func * jacobian[ix]
      self.g[ self.refinery_target.parameter_managers[ix].param_id ] = flex.sum(
         self.refinery_target.reflections[ "intensity_weight" ] * d)
      self.c[ self.refinery_target.parameter_managers[ix].param_id ] = flex.sum(
         self.refinery_target.reflections[ "intensity_weight" ] * d * d)
    print("rms %10.3f"%math.sqrt(flex.sum(self.refinery_target.reflections[ "intensity_weight" ]*self.func*self.func)/
                                 flex.sum(self.refinery_target.reflections[ "intensity_weight" ])), end=' ')#, file=self.out)
    print(list(self.x[-4:]))#,file=self.out)

    sel = self.c != 0
    diag = flex.double(len(self.c), 1)
    diag.set_selected(sel,  1 / self.c.select(sel))

    return self.f, self.g, diag


from scitbx.lstbx import normal_eqns, normal_eqns_solving

class per_frame_helper(normal_eqns.non_linear_ls, normal_eqns.non_linear_ls_mixin):
  def __init__(self,current_x=None, refinery_target=None, out=None,):
    self.outliers = False
    self.refinery_target = refinery_target
    self.out = out
    super(per_frame_helper, self).__init__(n_parameters=current_x.size())
    self.n = current_x.size()
    self.x_0 = current_x.deep_copy()
    self.restart()

  def restart(self):
    self.x = self.x_0.deep_copy()
    self.old_x = None

  def step_forward(self):
    self.old_x = self.x.deep_copy()
    self.x += self.step()

  def step_backward(self):
    assert self.old_x is not None
    self.x, self.old_x = self.old_x, None

  def parameter_vector_norm(self):
    return self.x.norm()

  def build_up(self, objective_only=False):
    P = Profiler("build up")
    # take current parameters self.x, translate into the data structure of the refinery_target
    self.refinery_target.take_current_values(self.x)
    residuals = self.refinery_target.fvec_callable()

    if not self.outliers:
      from scitbx.math import five_number_summary
      sel = self.refinery_target.reflections[ "intensity_weight" ] != 0
      print ("Initial count rejected", sel.count(False), "out of", len(sel))
      min_x, q1_x, med_x, q3_x, max_x = five_number_summary(residuals.select(sel))
      print ("Five number summary: min %.1f, q1 %.1f, med %.1f, q3 %.1f, max %.1f"%(min_x, q1_x, med_x, q3_x, max_x))
      iqr_x = q3_x - q1_x
      cut_x = 30 * iqr_x
      sel.set_selected(residuals > q3_x + cut_x, False)
      sel.set_selected(residuals < q1_x - cut_x, False)
      self.refinery_target.reflections[ "intensity_weight" ].set_selected(~sel, 0)
      print ("Rejecting", sel.count(False), "out of", len(sel))
      self.outliers = True

    self.reset()
    if objective_only:
      self.add_residuals(residuals, weights=self.refinery_target.reflections[ "intensity_weight" ])
    else:
      from scitbx import sparse
      dresiduals_dp = self.refinery_target.jacobian_callable()
      jacobian = sparse.matrix(len(self.refinery_target.reflections), self.n_parameters)
      for ix in range(len(self.refinery_target.parameter_managers)):
        gradient = dresiduals_dp[ix]
        sparse_column = sparse.matrix_column(len(gradient))
        sparse_column.set_selected(gradient!=0, gradient)
        sparse_block = sparse.matrix(len(gradient), 1)
        sparse_block[:,0] = sparse_column
        jacobian.assign_block(sparse_block, 0, self.refinery_target.parameter_managers[ix].param_id)

      self.add_equations(residuals, jacobian, weights=self.refinery_target.reflections[ "intensity_weight" ])
    print("rms %10.3f"%math.sqrt(flex.sum(self.refinery_target.reflections[ "intensity_weight" ]*residuals*residuals)/
                                 flex.sum(self.refinery_target.reflections[ "intensity_weight" ])), end=' ')#, file=self.out)
    logstr = ""
    for manager_id, manager in enumerate(self.refinery_target.parameter_managers[-4:]):
      manager_id += len(self.refinery_target.parameter_managers) - 4
      if manager.type == 'deff':
        logstr += manager.type + " %e "%(2/math.exp(self.x[manager_id]))
      elif manager.type == 'eta':
        logstr += manager.type + " %e "%math.exp(self.x[manager_id])
      else:
        logstr += manager.type + " % 14.10f "%self.x[manager_id]
    print(logstr, end = ' ')

class levmar_minimizer_ratchet(object):
  def __init__(self, current_x = None, refinery_target = None, out = None):
    self.helper = per_frame_helper(current_x = current_x,
                                   refinery_target = refinery_target,
                                   out = out)
    self.iterations = normal_eqns_solving.levenberg_marquardt_iterations(
      non_linear_ls = self.helper,
      track_all=True,
      gradient_threshold=1e-04,
      step_threshold=1e-04,
      tau=1e-03,
      n_max_iterations=200)

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
