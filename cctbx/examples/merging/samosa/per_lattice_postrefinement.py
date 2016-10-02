from __future__ import division
import sys
import math
from scitbx import matrix
from cctbx import miller
from dials.array_family import flex
from scitbx.math.tests.tst_weighted_correlation import simple_weighted_correlation

def legacy_cxi_merge_postrefinement(measurements_orig, params, i_model, miller_set, result):
    out = sys.stdout
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

    if True:
      from scitbx import lbfgs
      from libtbx import adopt_init_args

      pair1 = flex.int([pair[1] for pair in matches.pairs()])
      pair0 = flex.int([pair[0] for pair in matches.pairs()])
      #raw_input("go:")
      # narrow things down to the set that matches, only
      observations = observations.customized_copy(
        indices = flex.miller_index([observations.indices()[p] for p in pair1]),
        data = flex.double([observations.data()[p] for p in pair1]),
        sigmas = flex.double([observations.sigmas()[p] for p in pair1]),
      )
      observations_original_index = observations_original_index.customized_copy(
        indices = flex.miller_index([observations_original_index.indices()[p] for p in pair1]),
        data = flex.double([observations_original_index.data()[p] for p in pair1]),
        sigmas = flex.double([observations_original_index.sigmas()[p] for p in pair1]),
      )

      #IOBSVEC = flex.double([observations.data()[p] for p in pair1])
      #MILLER = flex.miller_index([observations_original_index.indices()[p] for p in pair1])
      IOBSVEC = observations.data()
      ICALCVEC = flex.double([i_model.data()[p] for p in pair0])
      MILLER = observations_original_index.indices()
      print >> out, "ZZZ",observations.size(), observations_original_index.size(), len(MILLER)
      ORI = result["current_orientation"][0]
      Astar = matrix.sqr(ORI.reciprocal_matrix())
      WAVE = result["wavelength"]
      BEAM = matrix.col((0.0,0.0,-1./WAVE))
      BFACTOR = 0.
      DSSQ = ORI.unit_cell().d_star_sq(MILLER)

      #calculation of correlation here
      class Empty:
        def __init__(CC): CC.n_obs=0; CC.n_rejected=0
      data = Empty()

      I_reference = flex.double([i_model.data()[pair[0]] for pair in matches.pairs()])
      I_observed = flex.double([observations.data()[pair[1]] for pair in matches.pairs()])
      use_weights = False # New facility for getting variance-weighted correlation

      if use_weights:
         #variance weighting
        I_weight = flex.double(
          [1./(observations.sigmas()[pair[1]])**2 for pair in matches.pairs()])
      else:
        I_weight = flex.double(len(observations.sigmas()), 1.)

      data.n_obs = len(matches.pairs())
      if params.include_negatives:
        data.n_rejected = 0
      else:
        non_positive = flex.bool(
          [observations.data()[pair[1]] <= 0 for pair in matches.pairs()])
        data.n_rejected = non_positive.count(True)
        I_weight.set_selected (non_positive, 0.)
      N = data.n_obs - data.n_rejected

      SWC = simple_weighted_correlation(I_weight, I_reference, I_observed)

      if params.postrefinement.algorithm=="rs":
        print("IN RS METHOD")
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

        class unpack(object):
         def __init__(YY,values):
          YY.reference = values # simply the flex double list of parameters
         def __getattr__(YY,item):
          if item=="thetax" : return YY.reference[3]
          if item=="thetay" : return YY.reference[4]
          if item=="G" :      return YY.reference[0]
          if item=="BFACTOR": return YY.reference[1]
          if item=="RS":      return YY.reference[2]
          return getattr(YY,item)

         def show(values):
          print >> out, "G: %10.7f"%values.G,
          print >> out, "B: %10.7f"%values.BFACTOR, \
                "RS: %10.7f"%values.RS, \
                "%7.3f deg %7.3f deg"%(
            180.*values.thetax/math.pi,180.*values.thetay/math.pi)

        def lorentz_callable(values):
          return get_partiality_array(values)

        def get_partiality_array(values):
          rs = values.RS
          Rh = get_Rh_array(values)
          rs_sq = rs*rs
          PB = rs_sq / ((2. * (Rh * Rh)) + rs_sq)
          return PB

      elif params.postrefinement.algorithm=="eta_deff":
        DVEC = ORI.unit_cell().d(MILLER)
        eta_init = 2. * result["ML_half_mosaicity_deg"][0] * math.pi/180.
        D_eff_init = 2.*result["ML_domain_size_ang"][0]
        current = flex.double([SWC.slope, BFACTOR, eta_init, 0., 0.,D_eff_init,])

        class unpack(object):
         def __init__(YY,values):
          YY.reference = values # simply the flex double list of parameters
         def __getattr__(YY,item):
          if item=="thetax" : return YY.reference[3]
          if item=="thetay" : return YY.reference[4]
          if item=="G" :      return YY.reference[0]
          if item=="BFACTOR": return YY.reference[1]
          if item=="ETA":      return YY.reference[2]
          if item=="DEFF":      return YY.reference[5]
          return getattr(YY,item)

         def show(values):
          print >> out, "%10.7f"%values.G,
          print >> out, "%10.7f"%values.BFACTOR, \
                "eta %10.7f"%values.ETA, \
                "Deff %10.2f"%values.DEFF, \
                "%7.3f deg %7.3f deg"%(
            180.*values.thetax/math.pi,180.*values.thetay/math.pi)

        def lorentz_callable(values):
          Rh = get_Rh_array(values)
          Rs = flex.double(len(MILLER),1./values.DEFF)+flex.double(len(MILLER),values.ETA/2.)/DVEC
          ratio = Rh / Rs
          ratio_abs = flex.abs(ratio)
          return ratio_abs

        def get_partiality_array(values):
          Rh = get_Rh_array(values)
          Rs = flex.double(len(MILLER),1./values.DEFF)+flex.double(len(MILLER),values.ETA/2.)/DVEC
          Rs_sq = Rs * Rs
          Rh_sq = Rh * Rh
          numerator = Rs_sq - Rh_sq
          denominator = values.DEFF * Rs * Rs_sq
          partiality = numerator / denominator
          return partiality

      def get_Rh_array(values):
        Rh = flex.double()
        eff_Astar = get_eff_Astar(values)
        for mill in MILLER:
          x = eff_Astar * matrix.col(mill)
          Svec = x + BEAM
          Rh.append(Svec.length() - (1./WAVE))
        return Rh

      def get_eff_Astar(values):
        thetax = values.thetax; thetay = values.thetay;
        effective_orientation = ORI.rotate_thru((1,0,0),thetax
           ).rotate_thru((0,1,0),thetay
           )
        return matrix.sqr(effective_orientation.reciprocal_matrix())

      def scaler_callable(values):
        PB = get_partiality_array(values)
        EXP = flex.exp(-2.*values.BFACTOR*DSSQ)
        terms = values.G * EXP * PB
        return terms

      def fvec_callable(values):
        PB = get_partiality_array(values)
        EXP = flex.exp(-2.*values.BFACTOR*DSSQ)
        terms = (values.G * EXP * PB * ICALCVEC - IOBSVEC)
        # Ideas for improvement
        #   straightforward to also include sigma weighting
        #   add extra terms representing rotational excursion: terms.concatenate(1.e7*Rh)
        return terms

      func = fvec_callable(unpack(current))
      functional = flex.sum(func*func)
      print >> out, "functional",functional

      class c_minimizer:

        def __init__(self, current_x=None,
                     min_iterations=0, max_calls=1000, max_drop_eps=1.e-5):
          adopt_init_args(self, locals())
          self.n = current_x.size()
          self.x = current_x
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
          values = unpack(self.x)
          assert -150. < values.BFACTOR < 150. # limits on the exponent, please
          self.func = fvec_callable(values)
          functional = flex.sum(self.func*self.func)
          self.f = functional
          DELTA = 1.E-7
          self.g = flex.double()
          for x in xrange(self.n):
            templist = list(self.x)
            templist[x]+=DELTA
            dvalues = flex.double(templist)

            dfunc = fvec_callable(unpack(dvalues))
            dfunctional = flex.sum(dfunc*dfunc)
            #calculate by finite_difference
            self.g.append( ( dfunctional-functional )/DELTA )
          self.g[2]=0.
          print >> out, "rms %10.3f"%math.sqrt(flex.mean(self.func*self.func)),
          values.show()
          return self.f, self.g

        def __del__(self):
          values = unpack(self.x)
          print >> out, "FINALMODEL",
          print >> out, "rms %10.3f"%math.sqrt(flex.mean(self.func*self.func)),
          values.show()

      try:
        MINI = c_minimizer( current_x = current )
      except AssertionError: # on exponential overflow
        return null_data(
               file_name=file_name, log_out=out.getvalue(), low_signal=True)
      scaler = scaler_callable(unpack(MINI.x))
      if params.postrefinement.algorithm=="rs":
        fat_selection = (lorentz_callable(unpack(MINI.x)) > 0.2)
      else:
        fat_selection = (lorentz_callable(unpack(MINI.x)) < 0.9)
      fat_count = fat_selection.count(True)

      #avoid empty database INSERT, if there are insufficient centrally-located Bragg spots:
      # in samosa, handle this at a higher level, but handle it somehow.
      #if fat_count < 3:
      #  return null_data(
      #         file_name=file_name, log_out=out.getvalue(), low_signal=True)
      print >> out, "On total %5d the fat selection is %5d"%(len(observations.indices()), fat_count)
      print >> out, "ZZZ",observations.size(), observations_original_index.size(), len(fat_selection), len(scaler)
      observations_original_index = observations_original_index.select(fat_selection)

      observations = observations.customized_copy(
        indices = observations.indices().select(fat_selection),
        data = (observations.data()/scaler).select(fat_selection),
        sigmas = (observations.sigmas()/scaler).select(fat_selection)
      )
      matches = miller.match_multi_indices(
        miller_indices_unique=miller_set.indices(),
        miller_indices=observations.indices())

      values = unpack(MINI.x)
      return get_eff_Astar(values), values.RS
