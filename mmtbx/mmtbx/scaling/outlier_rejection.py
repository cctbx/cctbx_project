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
import sys, os
import math
import string
from cStringIO import StringIO


class outlier_manager(object):
  def __init__(self,
               miller_obs,
               out=None):
    self.out=out
    if self.out is None:
      self.out=sys.stdout
    # the original miller array
    self.miller_obs = miller_obs

    if self.miller_obs.observation_type() is None:
      raise Sorry("Unknown observation type")

    # we make a working copy of the above miller array
    self.work_obs = self.miller_obs.deep_copy().set_observation_type(
      self.miller_obs )

    if not self.work_obs.is_xray_intensity_array():
      self.work_obs = self.work_obs.f_as_f_sq()

    if not self.miller_obs.is_xray_amplitude_array():
      self.miller_obs = self.miller_obs.f_sq_as_f()

    #-----------------------
    # These calculations are needed for wilson based outlier rejection
    #
    # Normalize the data
    normalizer = absolute_scaling.kernel_normalisation(
      self.work_obs, auto_kernel=True )
    self.norma_work = self.work_obs.customized_copy(
      data=normalizer.normalised_miller.data()/
           normalizer.normalised_miller.epsilons().data().as_double() )
    # split things into centric and acentric sets please
    self.centric_work = self.norma_work.select_centric().set_observation_type(
      self.norma_work)
    self.acentric_work = self.norma_work.select_acentric().set_observation_type(
      self.norma_work)

  def basic_wilson_outliers(self,
                            p_basic_wilson=1E-6,
                            return_data=False):
    p_acentric_single = 1.0 - (1.0 - flex.exp(-self.acentric_work.data() ))
    p_centric_single = 1.0 - erf(flex.sqrt(self.centric_work.data()/2.0) )

    acentric_selection = flex.bool( p_acentric_single > p_basic_wilson )
    centric_selection = flex.bool( p_centric_single > p_basic_wilson )

    # combine all in a single miller array
    all_flags = self.work_obs.customized_copy(
      indices = self.acentric_work.indices().concatenate(
                self.centric_work.indices() ),
      data = acentric_selection.concatenate( centric_selection )
      )
    all_p_values = self.work_obs.customized_copy(
      indices = self.acentric_work.indices().concatenate(
                self.centric_work.indices() ),
      data = p_acentric_single.concatenate( p_centric_single )
      )

    # get the order right
    all_flags = all_flags.common_set(self.miller_obs)
    all_p_values = all_p_values.common_set(self.miller_obs)

    # prepare a table with results please
    log_string = """
Outlier rejection based on basic Wilson statistics.
--------------------------------------------------

See Read, Acta Cryst. (1999). D55, 1759-1764. for details.
Reflections whose normalized intensity have an associated p-value
lower than %s are flagged as possible outliers.
    """%(p_basic_wilson)


    log_string = self.make_log_wilson(log_string, all_flags ,all_p_values )
    print >> self.out
    print >> self.out, log_string
    print >> self.out

    if not return_data:
      return all_flags
    else:
      return self.miller_obs.select( all_flags.data() )


  def extreme_wilson_outliers(self,
                              p_extreme_wilson=1e-1,
                              return_data=False):

    n_acentric = self.acentric_work.data().size()
    n_centric = self.centric_work.data().size()

    extreme_acentric = 1.0 -  \
       flex.pow(1.0 - flex.exp(-self.acentric_work.data() ),float(n_acentric))
    extreme_centric = 1.0 - \
       flex.pow(erf(flex.sqrt(self.centric_work.data()/2.0) ),float(n_centric))

    acentric_selection = flex.bool(extreme_acentric > p_extreme_wilson)
    centric_selection  = flex.bool(extreme_centric > p_extreme_wilson)
    all_flags = self.work_obs.customized_copy(
      indices = self.acentric_work.indices().concatenate(
                self.centric_work.indices() ),
      data    = acentric_selection.concatenate( centric_selection )
    )
    all_p_values = self.work_obs.customized_copy(
      indices = self.acentric_work.indices().concatenate(
                self.centric_work.indices() ),
      data = extreme_acentric.concatenate( extreme_centric )
      )
    all_flags = all_flags.common_set(self.miller_obs)
    all_p_values = all_p_values.common_set(self.miller_obs)


    log_string = """
Outlier rejection based on extreme value Wilson statistics.
-----------------------------------------------------------

Reflections whose normalized intensity have an associated p-value
lower than %s are flagged as possible outliers.
The p-value is obtained using extreme value distributions of the
Wilson distribution.
    """%(p_extreme_wilson)

    log_string = self.make_log_wilson(log_string, all_flags ,all_p_values )

    print >> self.out
    print >> self.out, log_string
    print >> self.out

    if not return_data:
      return all_flags
    else:
      return self.miller_obs.select( all_flags.data() )

  def beamstop_shadow_outliers(self,
                               level=0.01,
                               d_min=10.0,
                               return_data=False):

    # just make sure that things don't get to weerd
    assert level <= 0.3
    z_lim_ac = -math.log(1.0-level)
    z_select_ac = flex.bool( self.norma_work.data() <
                             z_lim_ac )
    # a first order approximation of the NZ of centric is
    # sqrt(2/pi)sqrt(z)
    z_lim_c = level*level*math.pi/2.0
    z_select_c = flex.bool( self.norma_work.data() <
                             z_lim_c )
    centric = self.norma_work.centric_flags().data()
    d_select = flex.bool( self.norma_work.d_spacings().data() >
                          d_min )

    # final selection: all conditions must be full filled
    # acentrics: centric:FALSE
    #            z_select_ac:TRUE
    #            d_select:TRUE
    #  centrics: centric:TRUE
    #            z_select_c:TRUE
    #            d_select:TRUE
    #
    #  final = acentric is True OR centric is True
    #  THIS SHOULD BE DONE WITH OPERATIONS ON FLEX ARRAYS
    #  AM GETTING VERY CONFUSED HOWEVER!
    #  THIS DIDN'T WORK :
    #  a_part = ~(~z_select_ac or ~d_select)
    #  a_part = ~( ~a_part or centric)
    #  c_part = ~(~z_select_c  or ~d_select)
    #  c_part = ~( ~c_part or ~centric)
    #
    #  final_selection = ~(a_part or c_part)
    tmp_final = []
    for zac, zc, ds, cf in zip(z_select_ac,
                               z_select_c,
                               d_select,
                               centric):
      if ds:
        if cf:
          if zc:
            tmp_final.append( False )
          else:
            tmp_final.append( True )
        if not cf:
          if zac:
            tmp_final.append( False )
          else:
            tmp_final.append( True )
      else:
        tmp_final.append( True )
    tmp_final = flex.bool( tmp_final )
    final_selection = self.norma_work.customized_copy(
      data = tmp_final
    )
    log_message = """
Possible outliers due to beamstop shadow
----------------------------------------

Reflection with normalized intensities lower than %4.3e (acentric)
or %4.3e (centric) and an associated d-spacing lower than %3.1f
are considered potential outliers.
The rationale is that these reflection could potentially be in
the shadow of the beamstop.
     """%(z_lim_ac, z_lim_c, d_min)

    final_selection = final_selection.map_to_asu()
    self.miller_obs = self.miller_obs.map_to_asu()
    final_selection = final_selection.common_set( self.miller_obs )
    assert final_selection.indices().all_eq( self.miller_obs.indices() )
    data = self.miller_obs.select( final_selection.data()
                                   ).set_observation_type( self.miller_obs )

    log_message = self.make_log_beam_stop( log_message,final_selection )
    print >> self.out, log_message


    if not return_data:
      return final_selection
    else:
      return( data )

  def model_based_outliers(self,
                           f_model,
                           alpha,
                           beta,
                           level=15,
                           return_data=False,
                           plot_out=None):

    self.miller_obs = self.miller_obs.map_to_asu()
    f_model = f_model.common_set(self.miller_obs).map_to_asu()
    alpha = alpha.common_set( f_model ).map_to_asu()
    beta = beta.common_set( f_model ).map_to_asu()
    # check if all is well
    assert self.miller_obs.indices().all_eq( f_model.indices() )
    assert self.miller_obs.indices().all_eq( alpha.indices() )
    assert self.miller_obs.indices().all_eq( beta.indices() )

    f_model_outlier_object = scaling.likelihood_ratio_outlier_test(
      self.miller_obs.data(),
      self.miller_obs.sigmas(),
      flex.abs( f_model.data() ),
      self.miller_obs.epsilons().data().as_double(),
      self.miller_obs.centric_flags().data(),
      alpha.data(),
      beta.data()
      )
    modes = f_model_outlier_object.posterior_mode()
    lik = f_model_outlier_object.log_likelihood()
    p_lik = f_model_outlier_object.posterior_mode_log_likelihood()
    s_der = f_model_outlier_object.posterior_mode_snd_der()

    ll_gain = p_lik-lik
    # The smallest vallue should be 0.
    # sometimes, due to numerical issues, it comes out
    # a wee bit negative. please repair that
    eps=1.0e-10
    zeros = flex.bool( ll_gain < eps )
    ll_gain = ll_gain.set_selected( zeros, eps )
    #use the ll_gain to computew p values
    p_values = ll_gain*2.0
    p_values = erf( flex.sqrt(p_values/2.0) )
    p_values = 1.0 - flex.pow( p_values, float(p_values.size()) )


    # select on likelihood
    # flags = f_model_outlier_object.flag_potential_outliers( level/2.0  )
    # or it mihgt be better to do it on p value
    flags = flex.bool(p_values > level )
    flags = self.miller_obs.customized_copy( data = flags )
    ll_gain = self.miller_obs.customized_copy( data = ll_gain )
    p_values = self.miller_obs.customized_copy( data = p_values )
    f_model = self.miller_obs.customized_copy( data = flex.abs(f_model.data() ) )

    s_der = f_model_outlier_object.posterior_mode_snd_der()
    s_der = flex.sqrt(-1.0/s_der)
    s_der = self.miller_obs.customized_copy( data = s_der )

    log_message = """

Model based outlier rejection.
------------------------------

Calculated amplitudes and estimated values of alpha and beta
are used to compute the log-likelihood of the observed amplitude.
The method is inspired by Read, Acta Cryst. (1999). D55, 1759-1764.
Outliers are rejected on the basis of the assumption that the log
likelihood differnce log[P(Fobs)]-log[P(Fmode)] is distributed
according to a Chi-square distribution
(see http://en.wikipedia.org/wiki/Likelihood-ratio_test ).
The outlier threshold of the p-value relates to the p-value of the
extreme value distribution of the chi-square distribution.

"""

    flags.map_to_asu()
    ll_gain.map_to_asu()
    p_values.map_to_asu()
    s_der.map_to_asu()

    assert flags.indices().all_eq( self.miller_obs.indices() )
    assert ll_gain.indices().all_eq( self.miller_obs.indices() )
    assert p_values.indices().all_eq( self.miller_obs.indices() )
    assert s_der.indices().all_eq( self.miller_obs.indices() )

    log_message = self.make_log_model( log_message,
                                       flags,
                                       ll_gain,
                                       p_values,
                                       f_model,
                                       alpha ,
                                       beta,
                                       plot_out)
    tmp_log=StringIO()
    print >> tmp_log, log_message
    # histogram of log likelihood gain values
    print >> tmp_log
    print >> tmp_log, "The histoghram of 2(LL-gain) values is shown below."
    print >> tmp_log, "  Note: 2(LL-gain) is approximately Chi-square distributed."
    print >> tmp_log
    print >> tmp_log, "  2(LL-gain)       Frequency"
    histo = flex.histogram( (p_lik-lik)*2.0, 15 )
    histo.show(f=tmp_log,format_cutoffs='%7.3f')

    print >> self.out, tmp_log.getvalue()

    if not return_data:
      return flags
    else:
      assert flags.indices().all_eq(  self.miller_obs.indices() )
      return self.miller_obs.select( flags.data() )

  def apply_scale_to_original_data(self,scale_factor,d_min=None):
    self.miller_obs = self.miller_obs.customized_copy(
      data = self.miller_obs.data()*scale_factor
      )
    if d_min is not None:
      self.miller_obs = self.miller_obs.resolution_filter(d_min = d_min)

  def make_log_wilson(self, log_message, flags, p_values):
    """ produces a 'nice' table of outliers and their reason for
    being an outlier using basic or extreme wilson statistics """

    header = ("Index", "E_obs", "Centric", "p-value" )
    flags = flags.common_set( self.norma_work)
    p_vals = p_values.common_set( self.norma_work )

    rogues = self.norma_work.select( ~flags.data() )
    p_vals = p_vals.select( ~flags.data() )

    rows = []
    table = "No outliers were found."

    for hkl,e,c,p in zip(rogues.indices(),
                         rogues.data(),
                         rogues.centric_flags().data(),
                         p_vals.data() ):
      this_row = [str(hkl), "%5.3f"%(math.sqrt(e)), str(c), "%5.3e"%(p) ]
      rows.append( this_row)
    if len(rows)>0:
      table = table_utils.format( [header]+rows,
                                  comments=None,
                                  has_header=True,
                                  separate_rows=False,
                                  prefix='| ',
                                  postfix=' |')
    final = log_message +"\n" + table
    return final

  def make_log_beam_stop(self,
                         log_message, flags):
    self.norma_work = self.norma_work.map_to_asu()
    self.miller_obs = self.miller_obs.map_to_asu()
    flags = flags.map_to_asu()




    data = self.miller_obs.select( ~flags.data() )
    evals = self.norma_work.select( ~flags.data() )



    header = ("Index", "d-spacing", "F_obs","E-value", "centric")
    table = "No outliers were found"
    rows = []
    if data.data().size() > 0:
      if data.data().size() < 500 :
        for hkl, d, fobs, e, c in zip(data.indices(),
                                      data.d_spacings().data(),
                                      data.data(),
                                      evals.data(),
                                      data.centric_flags().data() ):
          this_row = [str(hkl),
                      "%4.2f"%(d),
                      "%6.1f"%(fobs),
                      "%4.2f"%( math.sqrt(e) ),
                      str(c) ]
          rows.append( this_row )

        table = table_utils.format( [header]+rows,
                                    comments=None,
                                    has_header=True,
                                    separate_rows=False,
                                    prefix='| ',
                                    postfix=' |')
      else:
        table = """Over 500 outliers have been found."""


    final = log_message + "\n" + table +"\n \n"
    return final


  def make_log_model(self,
                     log_message, flags,
                     ll_gain, p_values,
                     f_model,
                     alpha, beta,plot_out=None):
    header = ("Index", "d-spacing", "F_obs", "F_model", "2(LL-gain)", "p-value", "alpha", "beta", "centric")
    table="No outliers were found"
    rows = []
    rogues = self.miller_obs.select( ~flags.data() )
    p_array = p_values.select( ~flags.data() )
    ll_array = ll_gain.select( ~flags.data() )
    fc_array = f_model.select( ~flags.data() )
    alpha_array = alpha.select( ~flags.data() )
    beta_array = beta.select( ~flags.data() )



    centric_flags = self.miller_obs.centric_flags().select( ~flags.data() )
    if rogues.indices().size() > 0:
      if rogues.indices().size() < 500 :
        sigmas = rogues.sigmas()
        if rogues.sigmas()==None:
          sigmas = rogues.d_spacings().data()*0+10.0

        for hkl, d, fo, fc, llg, p, a, b, c, s, e in zip(
          rogues.indices(),
          rogues.d_spacings().data(),
          rogues.data(),
          fc_array.data(),
          ll_array.data(),
          p_array.data(),
          alpha_array.data(),
          beta_array.data(),
          centric_flags.data(),
          sigmas,
          rogues.epsilons().data().as_double() ):

          this_row = [str(hkl),
                      "%4.2f"%(d),
                      "%6.1f"%(fo),
                      "%6.1f"%(fc),
                      "%5.2f"%(2*llg),
                      "%5.3e"%(p),
                      "%4.3f"%(a),
                      "%5.3e"%(b),
                      str(c) ]
          rows.append( this_row )

          if plot_out is not None:
            outlier_plots.plotit(
              fobs=fo,
              sigma=s,
              fcalc=fc,
              alpha=a,
              beta=b,
              epsilon=e,
              centric=c,
              out=plot_out,
              plot_title=str(hkl) )



        table = table_utils.format( [header]+rows,
                                    comments=None,
                                    has_header=True,
                                    separate_rows=False,
                                    prefix='| ',
                                    postfix=' |')

      else:
        table = "More then 500 outliers were found. This is very suspicious. Check data or limits."

    final = log_message + "\n" + table
    return final
