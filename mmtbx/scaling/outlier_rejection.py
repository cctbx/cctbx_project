from __future__ import division
from cctbx import miller
from mmtbx import scaling
from mmtbx.scaling import sigmaa_estimation
from cctbx.array_family import flex
from libtbx.utils import Sorry
import mmtbx.scaling
from mmtbx.scaling import absolute_scaling, outlier_plots
from scitbx.math import erf
from libtbx import table_utils
from libtbx.utils import null_out
import sys
import math
from libtbx.str_utils import StringIO

class outlier_manager(object):
  def __init__(self,
               miller_obs,
               r_free_flags,
               out=None):
    self.out=out
    if self.out is None:
      self.out=sys.stdout
    if out == "silent":
      self.out = null_out()

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

    self.r_free_flags = r_free_flags

    #-----------------------
    # These calculations are needed for wilson based outlier rejection
    #
    # Normalize the data
    normalizer = absolute_scaling.kernel_normalisation(
      self.work_obs, auto_kernel=True )
    self.norma_work = self.work_obs.customized_copy(
      data=normalizer.normalised_miller.data()/
           normalizer.normalised_miller.epsilons().data().as_double() )
    assert ( flex.min(self.norma_work.data()) >= 0 )
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
                           level=.01,
                           return_data=False,
                           plot_out=None):

    assert  self.r_free_flags is not None
    if(self.r_free_flags.data().count(True)==0):
      self.r_free_flags = self.r_free_flags.array(
        data = ~self.r_free_flags.data())
    sigmaa_estimator = sigmaa_estimation.sigmaa_estimator(
      miller_obs   = self.miller_obs,
      miller_calc  = f_model,
      r_free_flags = self.r_free_flags,
      kernel_width_free_reflections = 200,
      n_sampling_points = 20,
      n_chebyshev_terms = 13 )

    sigmaa_estimator.show(out=self.out)
    sigmaa = sigmaa_estimator.sigmaa()
    obs_norm = abs(sigmaa_estimator.normalized_obs)
    calc_norm = sigmaa_estimator.normalized_calc

    f_model_outlier_object = scaling.likelihood_ratio_outlier_test(
      f_obs=obs_norm.data(),
      sigma_obs=None,
      f_calc=calc_norm.data(),
      # the data is prenormalized, all epsies are unity
      epsilon=flex.double(calc_norm.data().size(), 1.0),
      centric=obs_norm.centric_flags().data(),
      alpha=sigmaa.data(),
      beta=1.0-sigmaa.data()*sigmaa.data()
      )
    modes = f_model_outlier_object.posterior_mode()
    lik = f_model_outlier_object.log_likelihood()
    p_lik = f_model_outlier_object.posterior_mode_log_likelihood()
    s_der = f_model_outlier_object.posterior_mode_snd_der()

    ll_gain = f_model_outlier_object.standardized_likelihood()

    # The smallest vallue should be 0.
    # sometimes, due to numerical issues, it comes out
    # a wee bit negative. please repair that
    eps=1.0e-10
    zeros = flex.bool( ll_gain < eps )
    p_values = ll_gain
    p_values = p_values.set_selected( zeros, eps )
    p_values = erf( flex.sqrt(p_values/2.0) )
    p_values = 1.0 - flex.pow( p_values, float(p_values.size()) )

    # select on p-values
    flags    = flex.bool(p_values > level )
    flags    = self.miller_obs.customized_copy( data = flags )
    ll_gain  = self.miller_obs.customized_copy( data = ll_gain )
    p_values = self.miller_obs.customized_copy( data = p_values )

    log_message = """

Model based outlier rejection.
------------------------------

Calculated amplitudes and estimated values of alpha and beta
are used to compute the log-likelihood of the observed amplitude.
The method is inspired by Read, Acta Cryst. (1999). D55, 1759-1764.
Outliers are rejected on the basis of the assumption that a scaled
log likelihood differnce 2(log[P(Fobs)]-log[P(Fmode)])/Q\" is distributed
according to a Chi-square distribution (Q\" is equal to the second
derivative of the log likelihood function of the mode of the
distribution).
The outlier threshold of the p-value relates to the p-value of the
extreme value distribution of the chi-square distribution.

"""

    flags.map_to_asu()
    ll_gain.map_to_asu()
    p_values.map_to_asu()

    assert flags.indices().all_eq( self.miller_obs.indices() )
    assert ll_gain.indices().all_eq( self.miller_obs.indices() )
    assert p_values.indices().all_eq( self.miller_obs.indices() )

    log_message = self.make_log_model( log_message,
                                       flags,
                                       ll_gain,
                                       p_values,
                                       obs_norm,
                                       calc_norm,
                                       sigmaa,
                                       plot_out)
    tmp_log=StringIO()
    print >> tmp_log, log_message
    # histogram of log likelihood gain values
    print >> tmp_log
    print >> tmp_log, "The histoghram of scaled (LL-gain) values is shown below."
    print >> tmp_log, "  Note: scaled (LL-gain) is approximately Chi-square distributed."
    print >> tmp_log
    print >> tmp_log, "  scaled(LL-gain)  Frequency"
    histo = flex.histogram( ll_gain.data(), 15 )
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
      ).set_observation_type( self.miller_obs )
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
      if e > 0:
        this_row = [str(hkl), "%5.3f"%(math.sqrt(e)), str(c), "%5.3e"%(p) ]
      else:
        this_row = [str(hkl), "%5.3f"%(0), str(c), " inf" ]
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
                     log_message,
                     flags,
                     ll_gain,
                     p_values,
                     e_obs,
                     e_calc,
                     sigmaa,
                     plot_out=None):
    header = ("Index", "d-spacing", "E_obs", "E_model", "Score", "p-value", "sigmaa", "centric")
    table="No outliers were found"
    rows = []
    rogues   = e_obs.select( ~flags.data() )
    p_array  = p_values.select( ~flags.data() )
    ll_array = ll_gain.select( ~flags.data() )
    ec_array = e_calc.select( ~flags.data() )
    sa_array = sigmaa.select( ~flags.data() )




    centric_flags = self.miller_obs.centric_flags().select( ~flags.data() )
    if rogues.indices().size() > 0:
      if rogues.indices().size() < 500 :
        sigmas = rogues.sigmas()
        if rogues.sigmas()==None:
          sigmas = rogues.d_spacings().data()*0+10.0

        for hkl, d, eo, ec, llg, p, sa, c, s, e in zip(
          rogues.indices(),
          rogues.d_spacings().data(),
          rogues.data(),
          ec_array.data(),
          ll_array.data(),
          p_array.data(),
          sa_array.data(),
          centric_flags.data(),
          sigmas,
          rogues.epsilons().data().as_double() ):

          this_row = [str(hkl),
                      "%4.2f"%(d),
                      "%6.3f"%(eo),
                      "%6.3f"%(ec),
                      "%5.2f"%(llg),
                      "%5.3e"%(p),
                      "%4.3f"%(sa),
                      str(c) ]
          rows.append( this_row )

          if plot_out is not None:
            outlier_plots.plotit(
              fobs=eo,
              sigma=s,
              fcalc=ec,
              alpha=sa,
              beta=1.0-sa*sa,
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
