from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from scitbx import matrix
from cctbx import miller
from dials.array_family import flex
from scitbx.math.tests.tst_weighted_correlation import simple_weighted_correlation
from libtbx import adopt_init_args
from xfel.cxi.postrefinement_legacy_rs import legacy_rs, rs_refinery, rs_parameterization, lbfgs_minimizer_base

def chosen_weights(observation_set, params):
    data = observation_set.data()
    sigmas = observation_set.sigmas()
    return {
      "unit": flex.double(len(data),1.),
      "variance": 1./(sigmas*sigmas),
      "gentle": flex.pow(flex.sqrt(flex.abs(data))/sigmas,2),
      "extreme": flex.pow(data/sigmas,2)
    } [ params.postrefinement.target_weighting ]

class updated_rs(legacy_rs):
  def __init__(self,measurements_orig, params, i_model, miller_set, result, out):
    measurements = measurements_orig.deep_copy()

    # Now manipulate the data to conform to unit cell, asu, and space group
    # of reference.  The resolution will be cut later.
    # Only works if there is NOT an indexing ambiguity!
    observations = measurements.customized_copy(
      anomalous_flag=not params.merge_anomalous,
      crystal_symmetry=miller_set.crystal_symmetry()
      ).map_to_asu()

    observations_original_index = measurements.customized_copy(
      anomalous_flag=not params.merge_anomalous,
      crystal_symmetry=miller_set.crystal_symmetry()
      )

    # Ensure that match_multi_indices() will return identical results
    # when a frame's observations are matched against the
    # pre-generated Miller set, self.miller_set, and the reference
    # data set, self.i_model.  The implication is that the same match
    # can be used to map Miller indices to array indices for intensity
    # accumulation, and for determination of the correlation
    # coefficient in the presence of a scaling reference.

    assert len(i_model.indices()) == len(miller_set.indices()) \
        and  (i_model.indices() ==
              miller_set.indices()).count(False) == 0

    matches = miller.match_multi_indices(
      miller_indices_unique=miller_set.indices(),
      miller_indices=observations.indices())

    pair1 = flex.int([pair[1] for pair in matches.pairs()])
    pair0 = flex.int([pair[0] for pair in matches.pairs()])
    # narrow things down to the set that matches, only
    observations_pair1_selected = observations.customized_copy(
      indices = flex.miller_index([observations.indices()[p] for p in pair1]),
      data = flex.double([observations.data()[p] for p in pair1]),
      sigmas = flex.double([observations.sigmas()[p] for p in pair1]),
    )
    observations_original_index_pair1_selected = observations_original_index.customized_copy(
      indices = flex.miller_index([observations_original_index.indices()[p] for p in pair1]),
      data = flex.double([observations_original_index.data()[p] for p in pair1]),
      sigmas = flex.double([observations_original_index.sigmas()[p] for p in pair1]),
    )
###################
    I_observed = observations_pair1_selected.data()
    chosen = chosen_weights(observations_pair1_selected, params)

    MILLER = observations_original_index_pair1_selected.indices()
    ORI = result["current_orientation"][0]
    Astar = matrix.sqr(ORI.reciprocal_matrix())
    WAVE = result["wavelength"]
    BEAM = matrix.col((0.0,0.0,-1./WAVE))
    BFACTOR = 0.

    #calculation of correlation here
    I_reference = flex.double([i_model.data()[pair[0]] for pair in matches.pairs()])
    I_invalid = flex.bool([i_model.sigmas()[pair[0]] < 0. for pair in matches.pairs()])
    use_weights = False # New facility for getting variance-weighted correlation

    if use_weights:
       #variance weighting
      I_weight = flex.double(
        [1./(observations_pair1_selected.sigmas()[pair[1]])**2 for pair in matches.pairs()])
    else:
      I_weight = flex.double(len(observations_pair1_selected.sigmas()), 1.)
    I_weight.set_selected(I_invalid,0.)
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
    """
    if params.include_negatives:
      SWC = simple_weighted_correlation(I_weight, I_reference, I_observed)
    else:
      non_positive = ( observations_pair1_selected.data() <= 0 )
      SWC = simple_weighted_correlation(I_weight.select(~non_positive),
            I_reference.select(~non_positive), I_observed.select(~non_positive))

    print("CORR: Old correlation is", SWC.corr, file=out)
    if params.postrefinement.algorithm=="rs2":
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
      refinery = rs2_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE,
        ICALCVEC = I_reference, IOBSVEC = I_observed, WEIGHTS = chosen)
      refinery.set_profile_shape(params.postrefinement.lineshape)

    func = refinery.fvec_callable(parameterization_class(current))
    functional = flex.sum(func*func)
    print("functional",functional, file=out)
    self.current = current; self.parameterization_class = parameterization_class
    self.refinery = refinery; self.out=out; self.params = params;
    self.miller_set = miller_set
    self.observations_pair1_selected = observations_pair1_selected;
    self.observations_original_index_pair1_selected = observations_original_index_pair1_selected
    self.i_model = i_model

  def run_plain(self):
    self.MINI = lbfgs_minimizer_derivatives( current_x = self.current,
        parameterization = self.parameterization_class, refinery = self.refinery,
        out = self.out )
    self.refined_mini = self.MINI

  def rs2_parameter_range_assertions(self,values):
    # New range assertions for refined variables
    assert 0 < values.G, "G-scale value out of range ( < 0 ) after rs2 refinement"
    assert -25 < values.BFACTOR and values.BFACTOR < 25, "B-factor value out of range ( |B|>25 ) after rs2 refinement"
    assert -0.5<180.*values.thetax/math.pi<0.5,"thetax value out of range ( |rotx|>.5 degrees ) after rs2 refinement"
    assert -0.5<180.*values.thetay/math.pi<0.5,"thetay value out of range ( |roty|>.5 degrees ) after rs2 refinement"

  def result_for_cxi_merge(self, file_name):
    values = self.get_parameter_values()
    self.rs2_parameter_range_assertions(values)
    scaler = self.refinery.scaler_callable(self.parameterization_class(self.MINI.x))

    partiality_array = self.refinery.get_partiality_array(values)
    p_scaler = flex.pow(partiality_array,
                        0.5*self.params.postrefinement.merge_partiality_exponent)

    fat_selection = (partiality_array > 0.2)
    fat_count = fat_selection.count(True)
    scaler_s = scaler.select(fat_selection)
    p_scaler_s = p_scaler.select(fat_selection)

    #avoid empty database INSERT, if insufficient centrally-located Bragg spots:
    # in samosa, handle this at a higher level, but handle it somehow.
    if fat_count < 3:
      raise ValueError("< 3 near-fulls after refinement")
    print("On total %5d the fat selection is %5d"%(
      len(self.observations_pair1_selected.indices()), fat_count), file=self.out)
    observations_original_index = \
      self.observations_original_index_pair1_selected.select(fat_selection)

    observations = self.observations_pair1_selected.customized_copy(
      indices = self.observations_pair1_selected.indices().select(fat_selection),
      data = (self.observations_pair1_selected.data().select(fat_selection)/scaler_s),
      sigmas = (self.observations_pair1_selected.sigmas().select(fat_selection)/(scaler_s * p_scaler_s))
    )
    matches = miller.match_multi_indices(
      miller_indices_unique=self.miller_set.indices(),
      miller_indices=observations.indices())

    I_weight = flex.double(len(observations.sigmas()), 1.)
    I_reference = flex.double([self.i_model.data()[pair[0]] for pair in matches.pairs()])
    I_invalid = flex.bool([self.i_model.sigmas()[pair[0]] < 0. for pair in matches.pairs()])
    I_weight.set_selected(I_invalid,0.)
    SWC = simple_weighted_correlation(I_weight, I_reference, observations.data())
    print("CORR: NEW correlation is", SWC.corr, file=self.out)
    print("ASTAR_FILE",file_name,tuple(self.refinery.get_eff_Astar(values)), file=self.out)
    self.final_corr = SWC.corr
    self.refined_mini = self.MINI
    #another range assertion
    assert self.final_corr > 0.1,"correlation coefficient out of range (<= 0.1) after rs2 refinement"

    return observations_original_index,observations,matches

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

  def __init__(self, current_x=None, parameterization=None, refinery=None, out=None,
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
    print("rms %10.3f"%math.sqrt(flex.sum(self.refinery.WEIGHTS*self.func*self.func)/
                                              flex.sum(self.refinery.WEIGHTS)), end=' ', file=self.out)
    values.show(self.out)
    return self.f, self.g
