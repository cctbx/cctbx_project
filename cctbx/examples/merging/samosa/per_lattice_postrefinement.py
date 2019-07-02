from __future__ import absolute_import, division, print_function
from six.moves import range
import math
from scitbx import matrix
from cctbx import miller
from dials.array_family import flex
from scitbx.math.tests.tst_weighted_correlation import simple_weighted_correlation
from libtbx import adopt_init_args, group_args

class legacy_cxi_merge_postrefinement(object):
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
    MILLER = observations_original_index_pair1_selected.indices()
    ORI = result["current_orientation"][0]
    Astar = matrix.sqr(ORI.reciprocal_matrix())
    WAVE = result["wavelength"]
    BEAM = matrix.col((0.0,0.0,-1./WAVE))
    BFACTOR = 0.

    #calculation of correlation here
    I_reference = flex.double([i_model.data()[pair[0]] for pair in matches.pairs()])
    use_weights = False # New facility for getting variance-weighted correlation

    if use_weights:
       #variance weighting
      I_weight = flex.double(
        [1./(observations_pair1_selected.sigmas()[pair[1]])**2 for pair in matches.pairs()])
    else:
      I_weight = flex.double(len(observations_pair1_selected.sigmas()), 1.)

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

    print("Old correlation is", SWC.corr, file=out)
    if params.postrefinement.algorithm=="rs":
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
      refinery = rs_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE,
        ICALCVEC = I_reference, IOBSVEC = I_observed)

    elif params.postrefinement.algorithm=="eta_deff":
      eta_init = 2. * result["ML_half_mosaicity_deg"][0] * math.pi/180.
      D_eff_init = 2.*result["ML_domain_size_ang"][0]
      current = flex.double([SWC.slope, BFACTOR, eta_init, 0., 0.,D_eff_init,])

      parameterization_class = eta_deff_parameterization
      refinery = eta_deff_refinery(ORI=ORI, MILLER=MILLER, BEAM=BEAM, WAVE=WAVE,
        ICALCVEC = I_reference, IOBSVEC = I_observed)

    func = refinery.fvec_callable(parameterization_class(current))
    functional = flex.sum(func*func)
    print("functional",functional, file=out)
    self.current = current; self.parameterization_class = parameterization_class
    self.refinery = refinery; self.out=out; self.params = params;
    self.miller_set = miller_set
    self.observations_pair1_selected = observations_pair1_selected;
    self.observations_original_index_pair1_selected = observations_original_index_pair1_selected

  def run_plain(self):
    self.MINI = lbfgs_minimizer_derivatives( current_x = self.current,
        parameterization = self.parameterization_class, refinery = self.refinery,
        out = self.out )

  def result_for_cxi_merge(self, file_name):
    scaler = self.refinery.scaler_callable(self.parameterization_class(self.MINI.x))
    if self.params.postrefinement.algorithm=="rs":
      fat_selection = (self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x)) > 0.2)
    else:
      fat_selection = (self.refinery.lorentz_callable(self.parameterization_class(self.MINI.x)) < 0.9)
    fat_count = fat_selection.count(True)

    #avoid empty database INSERT, if insufficient centrally-located Bragg spots:
    # in samosa, handle this at a higher level, but handle it somehow.
    if fat_count < 3:
      raise ValueError
    print("On total %5d the fat selection is %5d"%(
      len(self.observations_pair1_selected.indices()), fat_count), file=self.out)
    observations_original_index = \
      self.observations_original_index_pair1_selected.select(fat_selection)

    observations = self.observations_pair1_selected.customized_copy(
      indices = self.observations_pair1_selected.indices().select(fat_selection),
      data = (self.observations_pair1_selected.data()/scaler).select(fat_selection),
      sigmas = (self.observations_pair1_selected.sigmas()/scaler).select(fat_selection)
    )
    matches = miller.match_multi_indices(
      miller_indices_unique=self.miller_set.indices(),
      miller_indices=observations.indices())
    return observations_original_index,observations,matches

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
      Rh = flex.double()
      eff_Astar = self.get_eff_Astar(values)
      for mill in self.MILLER:
        x = eff_Astar * matrix.col(mill)
        Svec = x + self.BEAM
        Rh.append(Svec.length() - (1./self.WAVE))
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
      """
      hard_sphere_partial =(rs_sq-(Rh*Rh))/rs_sq
      immersion = Rh/rs
      from matplotlib import pyplot as plt
      plt.plot (immersion, PB, "r.")
      plt.plot (immersion, hard_sphere_partial,"b.")
      plt.plot ([-1.0,1.0],[0,0],"k-")
      plt.show()
      """
      return PB

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
      dPB_dRh = -PB * 4. * Rh / (2. * Rh * Rh + rs_sq)
      dPB_dthetax = dPB_dRh * dRh_dthetax
      dPB_dthetay = dPB_dRh * dRh_dthetay
      Px_terms = P_terms * dPB_dthetax; Py_terms = P_terms * dPB_dthetay

      return [G_terms,B_terms,0,Px_terms,Py_terms]

class nave1_refinery(rs_refinery):

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
      dPB_dRh = -PB * 4. * Rh / denomin
      dPB_dthetax = dPB_dRh * dRh_dthetax
      dPB_dthetay = dPB_dRh * dRh_dthetay
      Px_terms = P_terms * dPB_dthetax; Py_terms = P_terms * dPB_dthetay

      dPB_drs = 4 * rs * Rh * Rh / (denomin * denomin)
      Prs_terms = P_terms * dPB_drs

      return [G_terms,B_terms,Prs_terms,Px_terms,Py_terms]


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
  def __getattr__(YY,item):
    raise NotImplementedError
  def show(values,out):
    raise NotImplementedError

class rs_parameterization(unpack_base):
  def __getattr__(YY,item):
    if item=="thetax" : return YY.reference[3]
    if item=="thetay" : return YY.reference[4]
    if item=="G" :      return YY.reference[0]
    if item=="BFACTOR": return YY.reference[1]
    if item=="RS":      return YY.reference[2]
    return getattr(YY,item)

  def show(YY, out):
    print("G: %10.7f"%YY.G, end=' ', file=out)
    print("B: %10.7f"%YY.BFACTOR, \
        "RS: %10.7f"%YY.RS, \
        "%7.3f deg %7.3f deg"%(
        180.*YY.thetax/math.pi,180.*YY.thetay/math.pi), file=out)

class eta_deff_parameterization(unpack_base):
  def __getattr__(YY,item):
    if item=="thetax" : return YY.reference[3]
    if item=="thetay" : return YY.reference[4]
    if item=="G" :      return YY.reference[0]
    if item=="BFACTOR": return YY.reference[1]
    if item=="ETA":      return YY.reference[2]
    if item=="DEFF":      return YY.reference[5]
    return getattr(YY,item)

  def show(YY, out):
    print("%10.7f"%YY.G, end=' ', file=out)
    print("%10.7f"%YY.BFACTOR, \
          "eta %10.7f"%YY.ETA, \
          "Deff %10.2f"%YY.DEFF, \
          "%7.3f deg %7.3f deg"%(
      180.*YY.thetax/math.pi,180.*YY.thetay/math.pi), file=out)

class lbfgs_minimizer_base:

  def __init__(self, current_x=None, parameterization=None, refinery=None, out=None,
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
    self.g[2]=0.
    print("rms %10.3f"%math.sqrt(flex.mean(self.func*self.func)), end=' ', file=self.out)
    values.show(self.out)
    return self.f, self.g

  def __del__(self):
    values = self.parameterization(self.x)
    print("FINALMODEL", end=' ', file=self.out)
    print("rms %10.3f"%math.sqrt(flex.mean(self.func*self.func)), end=' ', file=self.out)
    values.show(self.out)

class lbfgs_minimizer_derivatives(lbfgs_minimizer_base):

  def compute_functional_and_gradients_test_code(self):
    values = self.parameterization(self.x)
    assert -150. < values.BFACTOR < 150. # limits on the exponent, please
    self.func = self.refinery.fvec_callable(values)
    functional = flex.sum(self.func*self.func)
    self.f = functional
    jacobian = self.refinery.jacobian_callable(values)
    self.gg_0 = flex.sum(2. * self.func * jacobian[0])
    self.gg_1 = flex.sum(2. * self.func * jacobian[1])
    self.gg_3 = flex.sum(2. * self.func * jacobian[3])
    self.gg_4 = flex.sum(2. * self.func * jacobian[4])
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
    self.g[2]=0.

    print("rms %10.3f"%math.sqrt(flex.mean(self.func*self.func)), end=' ', file=self.out)
    values.show(self.out)
    print("derivatives--> %15.5f    %15.5f    %9.7f   %5.2f   %5.2f"%tuple(self.g), file=self.out)
    print("  analytical-> %15.5f    %15.5f                %5.2f   %5.2f"%(
      self.gg_0,self.gg_1, self.gg_3,self.gg_4), file=self.out)
    self.g[0]=self.gg_0
    self.g[1]=self.gg_1
    self.g[3]=self.gg_3
    self.g[4]=self.gg_4
    return self.f, self.g
  def compute_functional_and_gradients(self):
    values = self.parameterization(self.x)
    assert -150. < values.BFACTOR < 150. # limits on the exponent, please
    self.func = self.refinery.fvec_callable(values)
    functional = flex.sum(self.func*self.func)
    self.f = functional
    jacobian = self.refinery.jacobian_callable(values)
    self.g = flex.double(self.n)
    for ix in range(self.n):
      self.g[ix] = flex.sum(2. * self.func * jacobian[ix])
    print("rms %10.3f"%math.sqrt(flex.mean(self.func*self.func)), end=' ', file=self.out)
    values.show(self.out)
    return self.f, self.g

