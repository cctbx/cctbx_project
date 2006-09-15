from __future__ import division
from cctbx import maptbx
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx import adptbx
from mmtbx import scaling
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.utils import Sorry, date_and_time, multi_out
from iotbx import reflection_file_reader
from iotbx import reflection_file_utils
from iotbx import crystal_symmetry_from_any
from iotbx import data_plots
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, outlier_plots
from scitbx.math import chebyshev_lsq
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit
from scitbx.math import erf
import libtbx.phil.command_line
from libtbx import table_utils
from scitbx.python_utils import easy_pickle
import scitbx.lbfgs
import sys, os
import math
import string
from cStringIO import StringIO
import mmtbx.f_model
from scitbx.math import chebyshev_lsq
from scitbx.math import chebyshev_polynome
from scitbx.math import chebyshev_lsq_fit



class sigmaa_point_estimator(object):
  def __init__(self,
               target_functor,
               h):
    self.functor = target_functor
    self.h = h
    self.f = None
    self.x = flex.double( [-1.0] )
    self.max = 0.99
    self.min = 0.01
    term_parameters = scitbx.lbfgs.termination_parameters(
      max_iterations = 100 )
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self,
                                      termination_params=term_parameters)

    self.sigmaa = self.min+(self.max-self.min)/(1.0+math.exp(-self.x[0]))

  def compute_functional_and_gradients(self):
    sigmaa = self.min+(self.max-self.min)/(1.0+math.exp(-self.x[0]))
    # chain rule bit for sigmoidal function
    dsdx = (self.max-self.min)*math.exp(-self.x[0])/(
      (1.0+math.exp(-self.x[0]))**2.0 )
    f = -self.functor.target(self.h,
                             sigmaa)
    g = -self.functor.dtarget(self.h,
                              sigmaa)*dsdx
    self.f = f
    return f, flex.double([g])



class sigmaa_estimator(object):
  def __init__(self,
               miller_obs,
               miller_calc,
               r_free_flags,
               width=None,
               number=100,
               auto_kernel=True,
               n_points=20,
               n_terms=13):
    self.n_terms=n_terms
    self.miller_obs = miller_obs.map_to_asu()
    self.miller_calc = abs(miller_calc.map_to_asu())
    self.r_free_flags = r_free_flags.map_to_asu()

    assert self.r_free_flags.indices().all_eq(
      self.miller_obs.indices() )

    self.miller_calc = self.miller_calc.common_set(
      self.miller_obs )
    assert self.r_free_flags.indices().all_eq(
      self.miller_calc.indices() )

    assert self.miller_obs.is_real_array()


    if self.miller_obs.is_xray_intensity_array():
      self.miller_obs = self.miller_obs.f_sq_as_f()
    assert self.miller_obs.is_xray_amplitude_array()
    self.miller_calc = self.miller_calc.set_observation_type(
      self.miller_obs)

    # get the aparameters that determine the kernel width
    self.width=width
    self.number=number
    self.auto_kernel=auto_kernel

    # get normalized data please
    self.normalized_obs_f = absolute_scaling.kernel_normalisation(
      self.miller_obs, auto_kernel=True)
    self.normalized_obs =self.normalized_obs_f.normalised_miller_dev_eps.f_sq_as_f()

    self.normalized_calc_f = absolute_scaling.kernel_normalisation(
      self.miller_calc, auto_kernel=True)
    self.normalized_calc =self.normalized_calc_f.normalised_miller_dev_eps.f_sq_as_f()

    # get the 'free data'
    self.free_norm_obs = self.normalized_obs.select( self.r_free_flags.data() )
    self.free_norm_calc= self.normalized_calc.select( self.r_free_flags.data() )


    # now we have to determin the width of the kernel we will use
    # use the same scheme as done for normalistion
    ## get the d_star_sq_array and sort it
    if self.auto_kernel:
      assert self.width is None
      assert self.number > 10 # be sensible please
      d_star_sq_hkl = self.free_norm_obs.d_star_sq().data()
      sort_permut = flex.sort_permutation(d_star_sq_hkl)

      if self.number > d_star_sq_hkl.size():
        self.number = d_star_sq_hkl.size()

      self.width = d_star_sq_hkl[sort_permut[number]]-flex.min( d_star_sq_hkl )
      assert self.width > 0

    if not self.auto_kernel:
      assert self.width is not None
      assert self.width > 0

    self.sigma_target_functor = scaling.sigmaa_estimator(
      e_obs     = self.free_norm_obs.data(),
      e_calc    = self.free_norm_calc.data(),
      centric   = self.free_norm_obs.centric_flags().data(),
      d_star_sq = self.free_norm_obs.d_star_sq().data() ,
      width=self.width)

    d_star_sq_overall = self.miller_obs.d_star_sq().data()
    self.min_h = flex.min( d_star_sq_overall )*0.99
    self.max_h = flex.max( d_star_sq_overall )*1.01
    self.h_array = flex.double( range(n_points) )*(
      self.max_h-self.min_h)/float(n_points-1.0)+self.min_h
    self.sigmaa_array = flex.double([])


    for h in self.h_array:
      stimator = sigmaa_point_estimator(self.sigma_target_functor,
                                        h)
      self.sigmaa_array.append( stimator.sigmaa )

    #now please fit a smooth  function to this bugger
    fit_lsq = chebyshev_lsq_fit.chebyshev_lsq_fit(
      self.n_terms,
      self.h_array,
      flex.log(self.sigmaa_array) )

    cheb_pol = chebyshev_polynome(
        self.n_terms,
        self.min_h,
        self.max_h,
        fit_lsq.coefs)

    self.sigmaa_array = flex.exp( cheb_pol.f( self.h_array ) )

    self.sigmaa = flex.exp( cheb_pol.f( d_star_sq_overall ) )
    self.alpha = self.sigmaa*flex.sqrt(
      self.normalized_calc_f.normalizer_for_miller_array/
      self.normalized_obs_f.normalizer_for_miller_array)
    self.beta = (1.0-self.sigmaa*self.sigmaa)*\
                self.normalized_calc_f.normalizer_for_miller_array
    # make them into miller arrays
    self.sigmaa = self.miller_obs.customized_copy(
      data = self.sigmaa )
    self.alpha = self.miller_obs.customized_copy(
      data = self.alpha )
    self.beta = self.miller_obs.customized_copy(
      data = self.beta )

    # now only FOM's need to be computed
    tmp_x = self.sigmaa.data()*self.normalized_calc.data()*self.normalized_obs.data()
    tmp_x = tmp_x / (1.0-self.sigmaa.data()*self.sigmaa.data())
    centric_fom = flex.tanh( tmp_x )
    acentric_fom = scitbx.math.bessel_i1_over_i0( 2.0*tmp_x )
    # we need to make sure centric and acentrics are not mixed up ...
    centrics = self.sigmaa.centric_flags().data()
    centric_fom  = centric_fom.set_selected( ~centrics, 0 )
    acentric_fom = acentric_fom.set_selected( centrics, 0 )
    final_fom =  centric_fom + acentric_fom
    self.fom = self.sigmaa.customized_copy(data=final_fom)



  def show(self, out=None):
    if out is None:
      out = sys.stdout
    print >> out
    print >> out, "SigmaA estimation summary"
    print >> out, "-------------------------"
    print >> out, "Kernel width  :  %6.2e"%(self.width)
    print >> out, "First N       :  %i"%(self.number)
    print >> out, "Auto kernel   : ", self.auto_kernel
    print >> out, "No. of points :  %i"%(self.h_array.size())
    print >> out, "No. of terms  :  %i"%(self.n_terms)
    print >> out
    print >> out, "1/d^2      d    sigmaA"
    for h,sa in zip( self.h_array, self.sigmaa_array):
      d = math.sqrt(1.0/h)
      print >> out, "%5.4f   %5.2f   %4.3f"%(h,d,sa)
    print >> out
    print >> out
