#ifndef MMTBX_SCALING_OUTLIER_H
#define MMTBX_SCALING_OUTLIER_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cctbx/miller/sym_equiv.h>
#include <scitbx/math/chebyshev.h>
#include <mmtbx/scaling/twinning.h>


namespace mmtbx { namespace scaling { namespace outlier{

  // Some model based outlier criteria are computed here
  template <typename FloatType>
  class likelihood_ratio_outlier_test{
  public:
    likelihood_ratio_outlier_test(scitbx::af::const_ref<FloatType> const& f_obs,
                                  scitbx::af::const_ref<FloatType> const& s_obs,
                                  scitbx::af::const_ref<FloatType> const& f_calc,
                                  scitbx::af::const_ref<FloatType> const& epsilon,
                                  scitbx::af::const_ref<bool> const& centric,
                                  scitbx::af::const_ref<FloatType> const& alpha,
                                  scitbx::af::const_ref<FloatType> const& beta ):
      log_ei0_(40000) // the exponentiated besseli0 function (  exp[-x]Io[x]  )
      {
        // make sure we have all equal sized arrays.
        SCITBX_ASSERT( f_obs.size() > 0 );
        SCITBX_ASSERT( f_calc.size() == f_obs.size() );
        if ( s_obs.size()!=0 ){
          SCITBX_ASSERT( s_obs.size() == f_obs.size() );
        }
        SCITBX_ASSERT( epsilon.size() == f_obs.size() );
        SCITBX_ASSERT( centric.size() == f_obs.size() );
        SCITBX_ASSERT( alpha.size() == f_obs.size() );
        SCITBX_ASSERT( beta.size() == f_obs.size() );

        // fill up the arrays please
        for (int ii=0; ii<f_obs.size();ii++){
          f_obs_.push_back( f_obs[ii] );
          if (s_obs.size() !=0){
            s_obs_.push_back( s_obs[ii] );
          } else {
            s_obs_.push_back( 0.0 );
          }
          f_calc_.push_back( f_calc[ii] );
          epsilon_.push_back( epsilon[ii] );
          centric_.push_back( centric[ii] );
          alpha_.push_back( alpha[ii]  );
          beta_.push_back( beta[ii] );
          // push back some zero's please
          f_obs_log_likelihood_.push_back(calc_log_likelihood(f_obs_[ii],ii) );
          f_obs_posterior_mode_.push_back(0.0);
          f_obs_posterior_mode_log_likelihood_.push_back(0.0);
          f_obs_posterior_mode_snd_der_.push_back(0.0);
          // get the stuff please
          newton_search(ii, 1.20, 1e-13, 1000);
        }
      }

      scitbx::af::shared<FloatType>
      log_likelihood(){
        return ( f_obs_log_likelihood_ );
      }

      scitbx::af::shared<FloatType>
      posterior_mode_log_likelihood(){
        return ( f_obs_posterior_mode_log_likelihood_ );
      }

      scitbx::af::shared<FloatType>
      posterior_mode(){
        return ( f_obs_posterior_mode_);
      }

      scitbx::af::shared<bool>
      flag_potential_outliers( FloatType level)
      {
        // the level of the test is the allowed difference
        // in log likelihood gain over the observed log likelihood and
        // that for the posterior mode (the ml estimate of Fobs given
        // alpha, beta and Fcalc).
        // When things were Gaussian, a level of 4.5 corresponds to
        // discarding differences larger than 3 sigma
        //
        scitbx::af::shared<bool> flags;
        FloatType delta;
        for (int ii=0;ii<f_calc_.size();ii++){
           delta = f_obs_posterior_mode_log_likelihood_[ii] - f_obs_log_likelihood_[ii];
           if (delta>level){
             flags.push_back( false );
           } else {
             flags.push_back( true );
           }
        }
        return( flags );
      }

  protected:
      // newton root finding
      // see
      // http://en.wikipedia.org/wiki/Newton%27s_method
      // for details of the algorithm
      // we use the updates
      // x_(n+1) = x_n - f'(x_n)/f"(x_n)
      // to find the maximum of a function
      // (a simple line search basically)
      void newton_search(int ii,
                         FloatType start_fraction,
                         FloatType eps,
                         int max_iter)
      {
        FloatType f_bar, start_f_bar;;
        FloatType f_bar_m1, f_bar_m2;
        FloatType funct, grad,delta,step;
        bool converged=false;
        int count=0;
        start_f_bar = f_calc_[ii];
        f_bar = f_calc_[ii]*start_fraction;
        f_bar_m1 = -1;
        f_bar_m2 = -1;
        while (!converged){
          // book keeping
          f_bar_m1 = f_bar;
          f_bar_m2 = f_bar_m1;
          // function and gradient
          funct = fst_der(f_bar, ii);
          grad  = snd_der(f_bar, ii);
          // newton step
          step = - funct/grad;

          // I don't want to overstep
          // allow a maximumstep size of
          if (step>0.1*start_f_bar){
            step=0.1*start_f_bar;
          }
          if (step < -0.1*start_f_bar){
            step=-0.1*start_f_bar;
          }

          f_bar = f_bar +step;
          // convergence tests
          delta = f_bar - f_bar_m2;
          if (delta<0){
            delta=-delta;
          }
          if (delta < eps){
            converged=true;
          }
          if (count >= max_iter){
            converged=true;
          }
          count++;
          //std::cout << count <<" " << f_bar << " " << funct << " " << grad << "  " << -funct/grad << std::endl;
        }

        f_obs_posterior_mode_log_likelihood_[ii]=calc_log_likelihood(f_bar,ii);
        f_obs_posterior_mode_[ii]=f_bar;
        f_obs_posterior_mode_snd_der_[ii]=grad;
      }

      inline FloatType fst_der(FloatType fo, int ii)
      {
        if (centric_[ii]){
          return( calc_fst_der_centric(fo,ii) );
        } else{
          return( calc_fst_der_acentric(fo,ii) );
        }
      }

      inline FloatType snd_der(FloatType fo, int ii)
      {
        if (centric_[ii]){
          return( calc_snd_der_centric(fo,ii) );
        } else{
          return( calc_snd_der_acentric(fo,ii) );
        }
      }

      // compute the first derivative of the loglikelihood function
      inline FloatType
      calc_fst_der_acentric( FloatType fo, int ii)
      {
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii] + s_obs_[ii]*s_obs_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = 2.0*alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = scitbx::math::bessel::i1_over_i0(x);
        result = (1.0/fo) - (2.0*fo/eb) +(2*alpha_[ii]*f_calc_[ii]/eb)*m;
        return (result);
      }

      // derivative of centrics
      inline FloatType
      calc_fst_der_centric( FloatType fo, int ii )
      {
        FloatType result;
        FloatType eb=2.0*epsilon_[ii]*beta_[ii] + s_obs_[ii]*s_obs_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = 2.0*alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = std::tanh(x)*2.0*alpha_[ii]*f_calc_[ii]/(eb);;
        result = -2.0*fo/eb + m;
        return (result);
      }


      // second derivatives
      inline FloatType
      calc_snd_der_acentric( FloatType fo, int ii)
      {
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii] + s_obs_[ii]*s_obs_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = 2.0*alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = scitbx::math::bessel::i1_over_i0(x);

        result = - (1.0/(fo*fo))
                 - (fo/eb)
                 - (2*alpha_[ii]/eb)*m
                 - (f_calc_[ii]*4.0*alpha_[ii]*alpha_[ii]/(eb*eb))*
                      (1.0-m/x-m*m);
        return (result);
      }

      inline FloatType calc_snd_der_centric( FloatType fo, int ii )
      {
        FloatType result;
        FloatType eb=2.0*epsilon_[ii]*beta_[ii] + s_obs_[ii]*s_obs_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = 2.0*alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = std::tanh( x );
        result = -(2.0/(eb)) -(2.0*alpha_[ii]*f_calc_[ii]/eb)*
          (2.0*alpha_[ii]*f_calc_[ii]/eb)*(1-m*m);
        return (result);
      }



      // likelihoods

      FloatType calc_log_likelihood(FloatType fo, int ii )
      {
        FloatType result=0;
        if (centric_[ii]){
          return( centric_log_likelihood(fo,ii) );
        } else {
          return( acentric_log_likelihood(fo,ii) );
        }
      }

      inline FloatType acentric_log_likelihood(FloatType fo, int ii){
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii] + s_obs_[ii]*s_obs_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x=2.0*alpha_[ii]*fo*f_calc_[ii]/eb;
        FloatType exparg = (fo - alpha_[ii]*f_calc_[ii]);
        exparg = exparg*exparg/eb;
        result = std::log(fo) - std::log(eb)  -exparg + log_ei0_.log_ei0(x);
        return (result);
      }


      inline FloatType centric_log_likelihood(FloatType fo, int ii){
        FloatType result;
        FloatType eb=2.0*epsilon_[ii]*beta_[ii] + s_obs_[ii]*s_obs_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x=2.0*alpha_[ii]*fo*f_calc_[ii]/eb;
        FloatType exparg = (fo*fo + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii])/eb;
        result = 0.5*std::log(2.0)-0.5*std::log(scitbx::constants::pi) - 0.5*std::log(0.5*eb)
          -exparg + std::log(std::cosh(x));
        return(result);
      }





      scitbx::af::shared<FloatType> f_obs_;
      scitbx::af::shared<FloatType> s_obs_;
      scitbx::af::shared<FloatType> f_calc_;
      scitbx::af::shared<FloatType> epsilon_;
      scitbx::af::shared<bool> centric_;
      scitbx::af::shared<FloatType> alpha_;
      scitbx::af::shared<FloatType> beta_;
      scitbx::af::shared<FloatType> f_obs_log_likelihood_;
      scitbx::af::shared<FloatType> f_obs_posterior_mode_;
      scitbx::af::shared<FloatType> f_obs_posterior_mode_log_likelihood_;
      scitbx::af::shared<FloatType> f_obs_posterior_mode_snd_der_;
      mmtbx::scaling::twinning::quick_log_ei0<FloatType> log_ei0_;

  };






}}}
#endif // MMTBX_SCALING_OUTLIER_H
