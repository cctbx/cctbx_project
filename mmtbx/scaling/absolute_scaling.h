//! Peter Zwart April 05, 2005
#ifndef MMTBX_SCALING_ABSOLUTE_SCALING_H
#define MMTBX_SCALING_ABSOLUTE_SCALING_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cctbx/miller/sym_equiv.h>
#include <mmtbx/scaling/scaling.h>

namespace mmtbx { namespace scaling{
namespace absolute_scaling{



  //! Negative of Wilson log likelihood for a single reflection
  /*!
   * Computes the negative log likelihood of a given
   * scale and wilson B value given a f_obs, sigma_f_obs
   * epsilon and the expected sum of squared form factors.
   *
   * FloatType const& f_obs        : F_obs
   * FloatType const& sigma_f_obs  : standard deviation of F_obs
   * FloatType const& epsilon      : epsilon
   * FloatType const& sig_sq       : sum of squared form factor
   * FloatType const& gamma        : a gamma term (see gamma_prot function)
   * bool centric                  : centric flag
   * FloatType const& p_scale      : (transformed) scale factor (scale=exp(-p_scale))
   * FloatType const& p_B_wilson   : wilson B value
   * bool transform                : whether or not scale and B should be transformed
   *                                 (default: True)
   */

  template <typename FloatType>
  FloatType
  wilson_single_nll(FloatType const& d_star_sq,
                    FloatType const& f_obs,
                    FloatType const& sigma_f_obs,
                    FloatType const& epsilon,
                    FloatType const& sig_sq,
                    FloatType const& gamma,
                    bool const& centric,
                    FloatType const& p_scale,
                    FloatType const& p_B_wilson,
                    bool const& transform=true)

  {

    SCITBX_ASSERT (f_obs>=0);
    SCITBX_ASSERT (sigma_f_obs>=0);

    typedef FloatType f_t;
    f_t B_wilson, scale=p_scale;
    if (transform){
      if (p_scale < -200.0){
        scale = -200.0;
      }
      if (p_scale > 200.0){
        scale=200.0;
      }
      scale = std::exp(-scale);
      B_wilson = p_B_wilson;
    } else {
      scale = scale;
      B_wilson = p_B_wilson;
    }

    f_t result = 0.0;
    f_t gamma_mult = 1.0+gamma;

    SCITBX_ASSERT (gamma_mult > 0); // this would be nonsense
    f_t k = std::max(1e-8,scale*std::exp(B_wilson*d_star_sq/4.0));
    f_t C = std::max(sig_sq*gamma_mult*epsilon + k*k*sigma_f_obs*sigma_f_obs,1e-8);
    if (!centric){
      result = - std::log(2.0)
               - std::log(k)
               - std::log(std::max(1e-12,f_obs))
               + std::log(C)
               + k*k*f_obs*f_obs/C;
      return(result);
    }
    result =  0.5*std::log(scitbx::constants::pi)
            + 0.5*std::log(C)
            + k*k*f_obs*f_obs/(2.0*C);
    return(result);
  }

  //! Computes the log likelihood for an array of reflections
  /*! See wilson_single_nll for details.
   *  Note that a resolution cut is made for safety.
   *  Ideally this should be done at a preprocessing level
   */
  template <typename FloatType>
  FloatType
  wilson_total_nll(scitbx::af::const_ref<FloatType> const& d_star_sq,
                   scitbx::af::const_ref<FloatType> const& f_obs,
                   scitbx::af::const_ref<FloatType> const& sigma_f_obs,
                   scitbx::af::const_ref<FloatType> const& epsilon,
                   scitbx::af::const_ref<FloatType> const& sig_sq,
                   scitbx::af::const_ref<FloatType> const& gamma,
                   scitbx::af::const_ref<bool> const& centric,
                   FloatType const& p_scale,
                   FloatType const& p_B_wilson,
                   bool transform=true)
  {

    SCITBX_ASSERT ( d_star_sq.size() == f_obs.size() );
    SCITBX_ASSERT ( d_star_sq.size() == sigma_f_obs.size() );
    SCITBX_ASSERT ( d_star_sq.size() == epsilon.size() );
    SCITBX_ASSERT ( d_star_sq.size() == sig_sq.size() );
    SCITBX_ASSERT ( d_star_sq.size() == gamma.size() );
    SCITBX_ASSERT ( d_star_sq.size() == centric.size() );


    typedef FloatType f_t;
    f_t result = 0.0;
    f_t tmp=0.0;
    unsigned n_data_points = d_star_sq.size();

    for (unsigned i=0; i<n_data_points; i++){
      if (d_star_sq[i] > d_star_sq_low_limit){
        if (d_star_sq[i] < d_star_sq_high_limit){
          tmp = wilson_single_nll(d_star_sq[i],
                                  f_obs[i],
                                  sigma_f_obs[i],
                                  epsilon[i],
                                  sig_sq[i],
                                  gamma[i],
                                  centric[i],
                                  p_scale,
                                  p_B_wilson,
                                  transform);
          result +=tmp;
        }
      }

    }
    return(result);

  }


  //! Gradients in ML isotropic wilson scaling for a single reflection
  /*!
   * Note that scale factors are allways transformed
   * in order to avoid problems with domain issues.
   * leave the domain for B from -inf to +inf
   * for the case that someone submits sharpened data
   */

  template <typename FloatType>
  scitbx::af::tiny<FloatType,2>
  wilson_single_nll_gradient(FloatType const& d_star_sq,
                             FloatType const& f_obs,
                             FloatType const& sigma_f_obs,
                             FloatType const& epsilon,
                             FloatType const& sig_sq,
                             FloatType const& gamma,
                             bool const& centric,
                             FloatType const& p_scale,
                             FloatType const& p_B_wilson)
  {

    typedef FloatType f_t;
    f_t EXP_ARG_MAX = 100.; // max value for exp argument to prevent from numerical issues
    f_t arg = -p_scale;
    if(arg > EXP_ARG_MAX) arg=EXP_ARG_MAX; // avoid overflow problem
    f_t kp=std::exp(arg);
    f_t Bp = p_B_wilson;
    f_t sigp = epsilon*sig_sq*(1.0+gamma);
    f_t sigd = sigma_f_obs*sigma_f_obs;
    f_t i_obs = f_obs*f_obs;
    arg = Bp*d_star_sq;
    if(arg > EXP_ARG_MAX) arg=EXP_ARG_MAX; // avoid overflow problem
    f_t exp_bp_ds_over_2 = std::exp(arg/2.0);
    f_t exp_bp_ds = std::exp(arg);
    f_t C = kp*kp*sigd*exp_bp_ds_over_2+sigp;
    CCTBX_ASSERT(C != 0.0); // avoid numerical issues
    f_t c_scale = 1./C;
    f_t c_scale_sq = c_scale * c_scale;
    f_t grad_scale = 0.0;
    f_t grad_B = 0.0;
    scitbx::af::tiny<double,2> gradients(0.0);
    if (centric){
      grad_scale = -exp_bp_ds*kp*kp*kp*i_obs*sigd*c_scale_sq
                 + exp_bp_ds_over_2*kp*i_obs*c_scale
                 + exp_bp_ds_over_2*kp*sigd*c_scale;
      grad_B = - exp_bp_ds*i_obs*d_star_sq*kp*kp*kp*kp*sigd/4.*c_scale_sq
               + exp_bp_ds_over_2*i_obs*d_star_sq*kp*kp/4.*c_scale
               + exp_bp_ds_over_2*d_star_sq*kp*kp*sigd/4.0*c_scale;
    } else {
      if(kp>1.e-9) {
        grad_scale = -1.0/kp
          - 2.0*exp_bp_ds*i_obs*kp*kp*kp*sigd*c_scale_sq
          + 2.0*exp_bp_ds_over_2*i_obs*kp*c_scale
          + 2.0*exp_bp_ds_over_2*kp*sigd*c_scale;
        grad_B = -d_star_sq/4.0
          - exp_bp_ds*i_obs*d_star_sq*kp*kp*kp*kp*sigd/2.0*c_scale_sq
          + exp_bp_ds_over_2*i_obs*d_star_sq*kp*kp/2.0*c_scale
          + exp_bp_ds_over_2*d_star_sq*kp*kp*sigd/2.0*c_scale;
      }
      else {
        grad_scale = 0.;
        grad_B = 0.;
      }
    }
    grad_scale = -grad_scale*kp;
    gradients[0] = grad_scale;
    gradients[1] = grad_B;
    return(gradients);
  }


  //!Compute the total gradient in ML Wilson scaling
  /*! Resolution cut is performed in this function
   *  for safety reasons, but should ideally be done at
   *  data preprocessing stage.
   */
  template <typename FloatType>
  scitbx::af::tiny<FloatType,2>
  wilson_total_nll_gradient(
    scitbx::af::const_ref<FloatType> const& d_star_sq,
    scitbx::af::const_ref<FloatType> const& f_obs,
    scitbx::af::const_ref<FloatType> const& sigma_f_obs,
    scitbx::af::const_ref<FloatType> const& epsilon,
    scitbx::af::const_ref<FloatType> const& sig_sq,
    scitbx::af::const_ref<FloatType> const& gamma,
    scitbx::af::const_ref<bool> const& centric,
    FloatType const& p_scale,
    FloatType const& p_B_wilson)
  {

    SCITBX_ASSERT ( d_star_sq.size() == f_obs.size() );
    SCITBX_ASSERT ( d_star_sq.size() == sigma_f_obs.size() );
    SCITBX_ASSERT ( d_star_sq.size() == epsilon.size() );
    SCITBX_ASSERT ( d_star_sq.size() == sig_sq.size() );
    SCITBX_ASSERT ( d_star_sq.size() == gamma.size() );
    SCITBX_ASSERT ( d_star_sq.size() == centric.size() );

    scitbx::af::tiny<FloatType,2> total_gradients(0,0);
    scitbx::af::tiny<FloatType,2> partial_gradients(0,0);

    unsigned n_data_points = d_star_sq.size();

    for (unsigned i=0;i<n_data_points;i++){
      if (d_star_sq[i] > d_star_sq_low_limit){
        if (d_star_sq[i] < d_star_sq_high_limit){
          partial_gradients =  wilson_single_nll_gradient(d_star_sq[i],
                                                          f_obs[i],
                                                          sigma_f_obs[i],
                                                          epsilon[i],
                                                          sig_sq[i],
                                                          gamma[i],
                                                          centric[i],
                                                          p_scale,
                                                          p_B_wilson);
          total_gradients[0]+=partial_gradients[0];
          total_gradients[1]+=partial_gradients[1];
        }
      }
    }

    return(total_gradients);
  }



  //! Normalisation of structure factors for given scale and B
  /*! Normalisation is done as follows:
   *  |E| = k*|F|/sigprot_sq
   *  Where
   *  k = exp[-p_scale]*exp[Bp*d_star_sq/4.0]
   *  and sigprot_sq is equal to the sum of squared structure factors
   *  Note that in this manner, one does not obtain <|E|**2>
   *  in each resolution bin. The <|E|**2> curve should follow
   *  the empirical gamma curve given in scaling.h
   *  This might be usefull for looking at sharpened maps
   *  or for a redetermination of an empirical gamma curve
   *
   *  Optionally (not the default), sigprot_sq is multiplied by 1+gamma
   *  These data can be used to perform various twinning tests if desired
   *
   *  scitbx::af::const_ref<FloatType> const& d_star_sq   : array of d_star_sq
   *  scitbx::af::const_ref<FloatType> const& f_obs       : array of f_obs's
   *  scitbx::af::const_ref<FloatType> const& epsilon     : array of epsilons
   *  scitbx::af::const_ref<FloatType> const& gamma       : array of gamma
   *  scitbx::af::const_ref<FloatType> const& sig_sq      : sig_prot_sq array
   *  scitbx::af::const_ref<bool> const& centric          : centric?
   *  FloatType const& p_scale                            : -log of scale
   *  FloatType const& p_B_wilson                         : isotropic B value
   *  bool wiggle                                         : if false:<|E|**1>=1
   *
   *  This function returns an const_Ref array with 'normalised'
   *  structure factor amplitudes.
   *
   *
   */

  template <typename FloatType>
  scitbx::af::shared<FloatType>
  ml_normalise(scitbx::af::const_ref<FloatType> const& d_star_sq,
               scitbx::af::const_ref<FloatType> const& f_obs,
               scitbx::af::const_ref<FloatType> const& epsilon,
               scitbx::af::const_ref<FloatType> const& sig_sq,
               scitbx::af::const_ref<FloatType> const& gamma,
               scitbx::af::const_ref<bool> const& centric,
               FloatType const& p_scale,
               FloatType const& p_B_wilson,
               bool const& wiggle=true)
  {

    SCITBX_ASSERT ( d_star_sq.size() == f_obs.size() );
    SCITBX_ASSERT ( d_star_sq.size() == epsilon.size() );
    SCITBX_ASSERT ( d_star_sq.size() == gamma.size() );
    SCITBX_ASSERT ( d_star_sq.size() == sig_sq.size() );
    SCITBX_ASSERT ( d_star_sq.size() == centric.size() );

    typedef FloatType f_t;
    f_t norma = 0;
    f_t scale = std::exp(-p_scale);
    f_t k;
    unsigned n_data_points = d_star_sq.size();
    scitbx::af::shared<FloatType> norma_sfa(n_data_points,0);

    for (unsigned i = 0; i<n_data_points; i++){
      // compute norma
      norma = sig_sq[i];
      if (!wiggle){
        norma = norma*(gamma[i]+1.0);
      }
      k = std::exp(p_B_wilson*d_star_sq[i]/4.0);
      norma_sfa[i] = f_obs[i]*k*scale/norma;
    }

    return(norma_sfa);
  }



  //! computing the scale factor (includes B-value correction)
  /*! See wilson_single_nll_aniso for details
   */

  template <typename FloatType>
  FloatType
  wilson_get_aniso_scale(cctbx::miller::index<> const& hkl,
                         FloatType const& p_scale,
                         FloatType const& V_star,
                         scitbx::sym_mat3<FloatType> const& u)
  {
    typedef FloatType f_t;

    f_t result = 0.0;

    result = hkl[0]*( hkl[0]*u[0]  + hkl[1]*u[3]  + hkl[2]*u[4] )
            +hkl[1]*( hkl[0]*u[3]  + hkl[1]*u[1]  + hkl[2]*u[5] )
            +hkl[2]*( hkl[0]*u[4]  + hkl[1]*u[5]  + hkl[2]*u[2] );

    result = result*scitbx::constants::pi
      *scitbx::constants::pi
      *2.0*V_star-p_scale;
      if(result>500) result=std::exp(500.);
      else result = std::exp(result);
    return (result);
  }



  //! Computing likelihood terms for anisotropic case
  /*! Ideally isotropic and anisotropic case should be combined
   *  in a single routine, and might happen after working out details.
   *
   * scitbx::af::const_ref<FloatType> const& H: Indices
   * FloatType const& f_obs                   : Observed amplitudes
   * FloatType const& sigma_f_obs             : Standard deviations
   * FloatType const& epsilon                 : Epsilons (multiplicity)
   * FloatType const& sig_sq                  : Squared of sum of form factors
   * FloatType const& gamma                   : Wilson wiggle
   * bool const& centric                      : Centric or acentric
   * FloatType const& p_scale                 : -logarithm of scale factor
   * cctbx::uctbx::unit_cell const& uc        : The unit cell
   * scitbx::af::const_ref<FloatType> U       : Scaled U tensor
   */

  template<typename FloatType>
  FloatType
  wilson_single_nll_aniso(cctbx::miller::index<> const& H,
                          FloatType const& f_obs,
                          FloatType const& sigma_f_obs,
                          FloatType const& epsilon,
                          FloatType const& sig_sq,
                          FloatType const& gamma,
                          bool const& centric,
                          FloatType const& p_scale,
                          cctbx::uctbx::unit_cell const& uc,
                          scitbx::sym_mat3<double> const& U)
  {
    SCITBX_ASSERT(H.size() == 3);
    SCITBX_ASSERT(U.size() == 6);

    typedef FloatType f_t;

    f_t result = 0.0;
    f_t k;
    f_t V_star_sq = pow( 1.0/(uc.volume()), 2.0/3.0);//RWGK's magic scalar

    k = wilson_get_aniso_scale(H,
                               p_scale,
                               V_star_sq,
                               U);
    f_t C=0;
    if(k<1.e+50 && sigma_f_obs<1.e+50) {
      C = epsilon*sig_sq*(1+gamma)+k*k*sigma_f_obs*sigma_f_obs;
    }
    if(k==0 || C==0 || C>1.e+50 || k>1.e+50) return 0;
    if (centric){
      result = 0.5*std::log(scitbx::constants::pi)
        + std::log(C)/2.0
        + k*k*f_obs*f_obs/(2.0*C);
      return(result);
    }
    result = -std::log(2.0)
      -std::log(k)
      -std::log(f_obs)
      +std::log(C)
      +k*k*f_obs*f_obs/C;
    return(result);
  }

  //!Computes the log likelihood for an array of reflections, anisotropic case
  /*!
   * See wilson_single_nll_aniso for details
   * Note that a resolution cut is made for safety.
   * Please take care of this at a preprocessing level.
   *
   */
  template<typename FloatType>
  FloatType
  wilson_total_nll_aniso(scitbx::af::const_ref<cctbx::miller::index<> >const&
                             hkl,
                         scitbx::af::const_ref<FloatType> const& f_obs,
                         scitbx::af::const_ref<FloatType> const& sigma_f_obs,
                         scitbx::af::const_ref<FloatType> const& epsilon,
                         scitbx::af::const_ref<FloatType> const& sig_sq,
                         scitbx::af::const_ref<FloatType> const& gamma,
                         scitbx::af::const_ref<bool> const& centric,
                         FloatType const& p_scale,
                         cctbx::uctbx::unit_cell const& uc,
                         scitbx::sym_mat3<FloatType> const& U)
  {

    SCITBX_ASSERT ( hkl.size() == f_obs.size() );
    SCITBX_ASSERT ( hkl.size() == sigma_f_obs.size() );
    SCITBX_ASSERT ( hkl.size() == epsilon.size() );
    SCITBX_ASSERT ( hkl.size() == sig_sq.size() );
    SCITBX_ASSERT ( hkl.size() == gamma.size() );
    SCITBX_ASSERT ( hkl.size() == centric.size() );

    typedef FloatType f_t;
    f_t result = 0.0;
    f_t tmp=0.0;
    f_t d_star_sq=0.0;

    unsigned n_data_points = hkl.size();

    for (unsigned i=0; i<n_data_points; i++){
      d_star_sq = uc.d_star_sq(hkl[i]);
      if (d_star_sq > d_star_sq_low_limit){
        if (d_star_sq < d_star_sq_high_limit){
          tmp = wilson_single_nll_aniso( hkl[i],
                                         f_obs[i],
                                         sigma_f_obs[i],
                                         epsilon[i],
                                         sig_sq[i],
                                         gamma[i],
                                         centric[i],
                                         p_scale,
                                         uc,
                                         U);
          result +=tmp;
        }

      }
    }

    return(result);
  }



  //! Computes the gradient for a single miller index
  template <typename FloatType>
  scitbx::af::shared<FloatType>
  wilson_single_nll_aniso_gradient(cctbx::miller::index<> const& hkl,
                                   FloatType const& f_obs,
                                   FloatType const& sigma_f_obs,
                                   FloatType const& epsilon,
                                   FloatType const& sig_sq,
                                   FloatType const& gamma,
                                   bool const& centric,
                                   FloatType const& p_scale,
                                   cctbx::uctbx::unit_cell const& uc,
                                   scitbx::sym_mat3<double> const& U)
    {
      SCITBX_ASSERT(hkl.size() == 3);
      SCITBX_ASSERT(U.size() == 6);

      typedef FloatType f_t;

      f_t k;
      f_t V_star = pow( 1.0/(uc.volume()), 2.0/3.0);//RWGK's magic scalar

      f_t grad_k=0.0;

      scitbx::af::shared<FloatType> gradient(7,0);

      k = wilson_get_aniso_scale(hkl,
                                 p_scale,
                                 V_star,
                                 U);
      f_t C = epsilon*sig_sq*(1+gamma)+k*k*sigma_f_obs*sigma_f_obs;
      if(k>1.e+50 || C>1.e+50 || C<1.e-50 || k<1.e-50) { // XXX Pavel: a quick fix to avoid overflow
        grad_k=0;                           // XXX Pavel: may want to look for the problem's root
      }
      else {
        if (centric){
          grad_k = -f_obs*f_obs*k*k*k*sigma_f_obs*sigma_f_obs/(C*C)
                   +f_obs*f_obs*k/C
                   +k*sigma_f_obs*sigma_f_obs/C;

        } else {
        grad_k = -(1.0)/(k)
                 -2.0*f_obs*f_obs*k*k*k*sigma_f_obs*sigma_f_obs/(C*C)
                 +2.0*f_obs*f_obs*k/C
                 +2.0*k*sigma_f_obs*sigma_f_obs/C;
        }
      }
      gradient[0] = -grad_k*k;
      f_t tmp = scitbx::constants::pi*scitbx::constants::pi*V_star;
      gradient[1]=2.0*tmp*hkl[0]*hkl[0]*k*grad_k;
      gradient[2]=2.0*tmp*hkl[1]*hkl[1]*k*grad_k;
      gradient[3]=2.0*tmp*hkl[2]*hkl[2]*k*grad_k;
      gradient[4]=4.0*tmp*hkl[0]*hkl[1]*k*grad_k;
      gradient[5]=4.0*tmp*hkl[0]*hkl[2]*k*grad_k;
      gradient[6]=4.0*tmp*hkl[1]*hkl[2]*k*grad_k;
      return(gradient);
    }



  //! This function computes the total gradient for the aniso. abs. scaling
  /*! Note that we loop over the complete orbit (i think that is what it is
   *  called) of the miller index h. Not reccomended as it is slowish.
   *  See wilson_total_nll_aniso_gradient for a better solution.
   */
  template<typename FloatType>
  scitbx::af::shared<FloatType>
  wilson_total_nll_aniso_gradient_orbit(
                         scitbx::af::const_ref<cctbx::miller::index<> >const&
                             hkl,
                         scitbx::af::const_ref<FloatType> const& f_obs,
                         scitbx::af::const_ref<FloatType> const& sigma_f_obs,
                         scitbx::af::const_ref<FloatType> const& epsilon,
                         scitbx::af::const_ref<FloatType> const& sig_sq,
                         scitbx::af::const_ref<FloatType> const& gamma,
                         scitbx::af::const_ref<bool> const& centric,
                         FloatType const& p_scale,
                         cctbx::uctbx::unit_cell const& uc,
                         cctbx::sgtbx::space_group const& space_group,
                         scitbx::sym_mat3<FloatType> const& U)
  {

    SCITBX_ASSERT ( hkl.size() == f_obs.size() );
    SCITBX_ASSERT ( hkl.size() == sigma_f_obs.size() );
    SCITBX_ASSERT ( hkl.size() == epsilon.size() );
    SCITBX_ASSERT ( hkl.size() == sig_sq.size() );
    SCITBX_ASSERT ( hkl.size() == gamma.size() );
    SCITBX_ASSERT ( hkl.size() == centric.size() );

    typedef FloatType f_t;

    scitbx::af::shared<FloatType> tmp(7,0);
    scitbx::af::shared<FloatType> result(7,0);

    f_t d_star_sq=0.25;

    unsigned n_data_points = hkl.size();


    f_t weight=0.0;
    cctbx::miller::index<> h_orbit;

    for (unsigned i=0; i<n_data_points; i++){

      d_star_sq = uc.d_star_sq(hkl[i]);

      if (d_star_sq > d_star_sq_low_limit){
        if (d_star_sq < d_star_sq_high_limit){

          cctbx::miller::sym_equiv_indices orbit(space_group, hkl[i]);
          weight = 1.0/orbit.indices().size();

          for (unsigned jj=0;jj<orbit.indices().size();jj++){
            h_orbit = orbit(jj).h();

            tmp = wilson_single_nll_aniso_gradient( h_orbit,
                                                    f_obs[i],
                                                    sigma_f_obs[i],
                                                    epsilon[i],
                                                    sig_sq[i],
                                                    gamma[i],
                                                    centric[i],
                                                    p_scale,
                                                    uc,
                                                    U);
            result[0] += tmp[0]*weight;
            result[1] += tmp[1]*weight;
            result[2] += tmp[2]*weight;
            result[3] += tmp[3]*weight;
            result[4] += tmp[4]*weight;
            result[5] += tmp[5]*weight;
            result[6] += tmp[6]*weight;
          }
        }
      }
    }
    return(result);
  }


  //! This function computes the total gradient for the aniso. abs. scaling
  /*! Here we do not loop over the orbit of hbut will do a post correction in
   *  python ( adp_constraints.independent_gradients() ).
   */
  template<typename FloatType>
  scitbx::af::shared<FloatType>
  wilson_total_nll_aniso_gradient(
                         scitbx::af::const_ref<cctbx::miller::index<> >const&
                             hkl,
                         scitbx::af::const_ref<FloatType> const& f_obs,
                         scitbx::af::const_ref<FloatType> const& sigma_f_obs,
                         scitbx::af::const_ref<FloatType> const& epsilon,
                         scitbx::af::const_ref<FloatType> const& sig_sq,
                         scitbx::af::const_ref<FloatType> const& gamma,
                         scitbx::af::const_ref<bool> const& centric,
                         FloatType const& p_scale,
                         cctbx::uctbx::unit_cell const& uc,
                         scitbx::sym_mat3<FloatType> const& U)
  {

    SCITBX_ASSERT ( hkl.size() == f_obs.size() );
    SCITBX_ASSERT ( hkl.size() == sigma_f_obs.size() );
    SCITBX_ASSERT ( hkl.size() == epsilon.size() );
    SCITBX_ASSERT ( hkl.size() == sig_sq.size() );
    SCITBX_ASSERT ( hkl.size() == gamma.size() );
    SCITBX_ASSERT ( hkl.size() == centric.size() );

    typedef FloatType f_t;

    scitbx::af::shared<FloatType> tmp(7,0);
    scitbx::af::shared<FloatType> result(7,0);

    f_t d_star_sq=0.25;

    unsigned n_data_points = hkl.size();




    for (unsigned i=0; i<n_data_points; i++){

      d_star_sq = uc.d_star_sq(hkl[i]);

      if (d_star_sq > d_star_sq_low_limit){
        if (d_star_sq < d_star_sq_high_limit){
          tmp = wilson_single_nll_aniso_gradient( hkl[i],
                                                  f_obs[i],
                                                  sigma_f_obs[i],
                                                  epsilon[i],
                                                  sig_sq[i],
                                                  gamma[i],
                                                  centric[i],
                                                  p_scale,
                                                  uc,
                                                  U);
          result[0] += tmp[0];
          result[1] += tmp[1];
          result[2] += tmp[2];
          result[3] += tmp[3];
          result[4] += tmp[4];
          result[5] += tmp[5];
          result[6] += tmp[6];
        }
      }
    }
    return(result);
  }


  //! Apply scale and anistropic B-value correction on the data
  /*! Use anisotropic tenmsor determined previously
   *  or something similar. Note that the RWGK scaling by V_star**(-2/3)
   *  of U_star is optional.
   *  It is howeve important to supply a U_star type anisotropic
   *  tensor. Providing U_cif, U_cart or B_cart will give meansingless results.
   *
   *  Note that this method does not provide means to get rid of the
   *  wilson wiggles, as the istopic variant does.
   *  It is better (and more convenient) to use the chebyshev functions
   *  for that purpose, see resolution_functions.h and related files.
   *
   */
  template <typename FloatType>
  scitbx::af::shared<FloatType>
  ml_normalise_aniso(scitbx::af::const_ref<cctbx::miller::index<> >const&
                       hkl,
                       scitbx::af::const_ref<FloatType> const& f_obs,
                       FloatType const& p_scale,
                       cctbx::uctbx::unit_cell const& uc,
                       scitbx::sym_mat3<FloatType> const& U,
                       bool const& volume_correction_rwgk=false)
  {

    SCITBX_ASSERT ( hkl.size() == f_obs.size() );

    scitbx::af::shared<FloatType> normalised(hkl.size(),0);
    FloatType k;
    FloatType V_star = 1.0;
    if (volume_correction_rwgk){
      V_star = pow( 1.0/(uc.volume()), 2.0/3.0);//RWGK's magic scalar
    }
    for (unsigned ii=0; ii<hkl.size(); ii++){
      k = wilson_get_aniso_scale(hkl[ii],
                                 p_scale,
                                 V_star,
                                 U);

      normalised[ii] = f_obs[ii]*k;
    }

    return (normalised);
  }


  template <typename FloatType>
  scitbx::af::shared<FloatType>
  kernel_normalisation(scitbx::af::const_ref<FloatType> const& d_star_sq_hkl,
                       scitbx::af::const_ref<FloatType> const& I_hkl,
                       scitbx::af::const_ref<FloatType> const& epsilon_hkl,
                       scitbx::af::const_ref<FloatType> const& d_star_sq_array,
                       FloatType const& kernel_width)
  {
    SCITBX_ASSERT(d_star_sq_hkl.size() == I_hkl.size() );
    SCITBX_ASSERT(d_star_sq_hkl.size() == epsilon_hkl.size() );

    scitbx::af::shared<FloatType> norma_I_array(d_star_sq_array.size(),0);
    scitbx::af::shared<FloatType> weights_array(d_star_sq_array.size(),0);
    FloatType x,dx,tmp_norm,result;
    FloatType eps=1E-8;

    // Use a simple kernel for binning purposes
    for (unsigned jj=0;jj<d_star_sq_hkl.size();jj++){
      x = d_star_sq_hkl[jj];
      for (unsigned ii=0 ;ii<d_star_sq_array.size();ii++){
        dx = x-d_star_sq_array[ii];
        dx = (dx*dx)/(2.0*kernel_width*kernel_width);
        result = std::exp(-dx);
        weights_array[ii] +=  result;
        norma_I_array[ii] +=  I_hkl[jj]*result/epsilon_hkl[jj];
      }
    }
    // Now we just have to 'normalise' via the weights we have obtained

    for (unsigned ii=0 ;ii<d_star_sq_array.size();ii++){
      tmp_norm = weights_array[ii];
      if (tmp_norm <= eps){
        tmp_norm=eps;
      }
      norma_I_array[ii]/=tmp_norm;
    }

    return(norma_I_array);
  }







}}}  // namespace mmtbx::scaling::absolute_scaling
#endif // MMTBX_SCALING_ABSOLUTE_SCALING_H
