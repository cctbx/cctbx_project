// FIXME this breaks in quite a few places on specific datasets - 1pcq in
// particular will fail several assertions.

#ifndef MMTBX_SCALING_OUTLIER_H
#define MMTBX_SCALING_OUTLIER_H

#include <scitbx/constants.h>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <cctbx/miller/sym_equiv.h>
#include <scitbx/math/chebyshev.h>
#include <mmtbx/scaling/twinning.h>
#include <scitbx/math/quadrature.h>
#include <scitbx/math/bessel.h>
#include <scitbx/line_search/more_thuente_1994.h>
#include <mmtbx/error.h>


namespace mmtbx { namespace scaling { namespace outlier{

  // Some model based outlier criteria are computed here
  // This function is tested having in mind that
  // alpha doesn't contain any scale factor.
  template <typename FloatType>
  class likelihood_ratio_outlier_test{
  public:
    likelihood_ratio_outlier_test(
          scitbx::af::const_ref<FloatType> const& f_obs,
          scitbx::af::const_ref<FloatType> const& s_obs,
          scitbx::af::const_ref<FloatType> const& f_calc,
          scitbx::af::const_ref<FloatType> const& epsilon,
          scitbx::af::const_ref<bool>      const& centric,
          scitbx::af::const_ref<FloatType> const& alpha,
          scitbx::af::const_ref<FloatType> const& beta )
      {
        FloatType min_beta=1.0e-5;
        FloatType min_alpha=1.0e-3;//
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

          if (alpha[ii]> min_alpha){
            alpha_.push_back( alpha[ii]  );
          } else {
             alpha_.push_back( min_alpha  );
          }
          if (beta[ii]> min_beta){
            beta_.push_back( beta[ii] );
          } else{
            beta_.push_back( min_beta );
          }
          // push back some zero's please
          f_obs_log_likelihood_.push_back(
            calc_log_likelihood(f_obs_[ii], ii) );
          f_obs_posterior_mode_.push_back(0.0);
          f_obs_posterior_mode_log_likelihood_.push_back(0.0);
          f_obs_posterior_mode_snd_der_.push_back(0.0);
          f_obs_fst_der_.push_back( fst_der(f_obs_[ii], ii)  );
          f_obs_snd_der_.push_back( snd_der(f_obs_[ii], ii)  );
          mean_fo_.push_back( compute_mean_( ii ) );
          std_fo_.push_back( compute_sigma_( ii )  );
          // get the stuff please
          newton_search(ii, 1.50,  1e-5, 500);
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

      scitbx::af::shared<FloatType>
      posterior_mode_snd_der(){
        return(f_obs_posterior_mode_snd_der_);
      }

      scitbx::af::shared<FloatType>
      f_obs_snd_der(){
        return(f_obs_snd_der_);
      }

      scitbx::af::shared<FloatType>
      f_obs_fst_der(){
        return(f_obs_fst_der_);
      }

      scitbx::af::shared<FloatType>
      mean_fobs(){
        return( mean_fo_  );
      }
      scitbx::af::shared<FloatType>
      std_fobs(){
        return( std_fo_ );
      }



      scitbx::af::shared<bool>
      flag_potential_outliers( FloatType level)
      {
        // the level of the test is the allowed difference
        // in log likelihood gain over the observed log likelihood and
        // that for the posterior mode (the ml estimate of Fobs given
        // alpha, beta and Fcalc.
        // When thnigs were Gaussian, a level of 4.5 corresponds to
        // discarding differences larger than 3 sigma
        //
        scitbx::af::shared<bool> flags;
        FloatType delta;
        for (int ii=0;ii<f_calc_.size();ii++){
           delta = f_obs_posterior_mode_log_likelihood_[ii] -
                   f_obs_log_likelihood_[ii];
           if ( (delta)>level){
             flags.push_back( false );
           } else {
             flags.push_back( true );
           }
        }
        return( flags );
      }


      scitbx::af::shared<FloatType>
      standardized_likelihood()
      {
        scitbx::af::shared<FloatType> result;
        FloatType delta;
        for (int ii=0;ii<f_obs_.size();ii++){
           delta = f_obs_posterior_mode_log_likelihood_[ii] -
                   f_obs_log_likelihood_[ii];
           result.push_back( 2.0*(delta) );
        }
        return( result );
      }




  protected:
      void newton_search(int ii,
                         FloatType start_fraction,
                         FloatType eps,
                         int max_iter)
      {
        // do a proper more-thuente line search
        FloatType eps_step=std_fo_[ii]/100.0, eps_delta=1e-9;
        scitbx::line_search::more_thuente_1994<FloatType>
          line_search_object;
        line_search_object.gtol=0.15;
        scitbx::af::shared<FloatType> xx, fp, sd;
        xx.push_back(0);
        scitbx::af::ref<FloatType> x = xx.ref();
        fp.push_back(0);
        sd.push_back(0);

        FloatType f, fpp, step_length, f_bar, f_bar_last=-10, delta;
        bool converged=false;
        int code;

        f_bar = mean_fo_[ii];
        int count=0;
        while (!converged){
          // restart the line search
          x[0]       =  f_bar;
          f_bar_last =  f_bar;
          f          = -calc_log_likelihood(x[0],ii);
          fp[0]      = -fst_der(x[0],ii);
          sd[0]      = -fp[0];
          fpp        = -snd_der(x[0],ii);
          //std::cout << "-------------" << std::endl;
          ///std::cout << count << " " << x[0] << " " << f << " " << f_calc_[ii]
          //        << " " << alpha_[ii] << " " << beta_[ii] << std::endl;
          //std::cout << fp[0] << std::endl;
          if  ( std::fabs(fp[0]) <= eps ){
            converged=true;
          }
          if (fpp<=eps_step){
            step_length = eps_step;
          } else {
            MMTBX_ASSERT(fpp != 0);
            step_length = 1.0/fpp;
          }
          if (!converged){
            code = line_search_object.start(x,
                                            f,
                                            fp.const_ref(),
                                            sd,
                                            step_length
                                            );
            if ( code == -1 ){
              while ( line_search_object.info_code==-1 ){
                f =-calc_log_likelihood( x[0], ii );
                fp[0]= fst_der( x[0],ii );
                code = line_search_object.next(x,
                                               f,
                                               fp.const_ref()
                                               );
              }
            } else {
              converged = true;
            }
          }

          f_bar = x[0];
          delta = std::fabs( f_bar - f_bar_last);
          if (delta < eps_delta){
            converged=true;
          }
          if (count>max_iter){
            converged=true;
          }
          count++;
        }

        FloatType grad;
        f_bar = x[0];
        if (f_bar < 0){
          f_bar = 0.0;
        }
        grad  = snd_der( f_bar, ii );
        FloatType tmp_rat = 2.0;
        if ( 1.0/std::sqrt(std::fabs(grad)) > tmp_rat*std_fo_[ii]){
          grad = std_fo_[ii];
          MMTBX_ASSERT(grad != 0);
          grad = -1.0/(grad*grad);
        }
        if ( grad > 0 ){
          grad = std_fo_[ii];
          MMTBX_ASSERT(grad != 0);
          grad = -1.0/(grad*grad);
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
        FloatType eb=epsilon_[ii]*beta_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        MMTBX_ASSERT((eb != 0) && (fo != 0));
        FloatType x = 2.0*alpha_[ii]*fo*f_calc_[ii]/(eb);

        FloatType m = scitbx::math::bessel::i1_over_i0(x);

        result = (1.0/fo) - (2.0*fo/eb) +(2.0*alpha_[ii]*f_calc_[ii]/eb)*m;
        return (result);
      }

      // derivative of centrics
      inline FloatType
      calc_fst_der_centric( FloatType fo, int ii )
      {
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii];
        MMTBX_ASSERT(eb != 0);
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = alpha_[ii]*fo*f_calc_[ii]/(eb);
        if (x<1e-13){
          x=1e-13;
        }
        FloatType m = std::tanh(x)*alpha_[ii]*f_calc_[ii]/(eb);;
        result = -fo/eb + m;
        return (result);
      }


      // second derivatives
      inline FloatType
      calc_snd_der_acentric( FloatType fo, int ii)
      {
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii];
        MMTBX_ASSERT(eb != 0);
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = 2.0*alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = scitbx::math::bessel::i1_over_i0(x);
        if (x<1e-13){
          x = 1e-13;
        }
        MMTBX_ASSERT((fo != 0) && (x != 0));
        result = - (1.0/(fo*fo))
                 - (2.0/eb)
                 + (f_calc_[ii]*4.0*alpha_[ii]*alpha_[ii]/(eb*eb))*
          (1.0-m/(x)-m*m);
        return (result);
      }

      inline FloatType calc_snd_der_centric( FloatType fo, int ii )
      {
        FloatType result=0;
        FloatType eb=epsilon_[ii]*beta_[ii];
        MMTBX_ASSERT(eb != 0);
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = std::tanh( x );
        result = -1.0/eb + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii]*(1.0-m*m)/(eb*eb);
        return (result);
      }



      // likelihoods

      FloatType calc_log_likelihood(FloatType fo, int ii )
      {
        if (centric_[ii]){
          return( centric_log_likelihood(fo,ii) );
        } else {
          return( acentric_log_likelihood(fo,ii) );
        }
      }

      inline FloatType acentric_log_likelihood(FloatType fo, int ii){
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        MMTBX_ASSERT(eb != 0);
        FloatType x=2.0*alpha_[ii]*fo*f_calc_[ii]/eb;
        FloatType exparg; // = (fo*fo + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii]);
        //exparg = exparg/eb;
        //result = std::log(2.0) + std::log(fo) - std::log(eb)  -exparg  + std::log( scitbx::math::bessel::i0(x) );
        //FloatType result2;
        exparg = fo -  alpha_[ii]*f_calc_[ii];
        exparg = exparg*exparg/eb;
        FloatType ei0x = scitbx::math::bessel::ei0(x);
        MMTBX_ASSERT((ei0x != 0) && (fo != 0));
        result = std::log(2.0) + std::log(fo) - std::log(eb)  -exparg + \
          std::log(ei0x);
        return (result);
      }


      inline FloatType centric_log_likelihood(FloatType fo, int ii){
        FloatType result;
        FloatType eb=epsilon_[ii]*beta_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        MMTBX_ASSERT(eb != 0);
        FloatType x=alpha_[ii]*fo*f_calc_[ii]/eb;
        FloatType exparg = (fo*fo + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii])/(2.0*eb);
        FloatType tmp;
        if (x>40){
           tmp = x*0.999921 - 0.65543;
        } else {
           tmp = std::log( std::cosh(x) );
        }
        result = 0.5*std::log(2.0)-0.5*std::log(scitbx::constants::pi) - 0.5*std::log(eb)
          -exparg + tmp;
        return(result);
      }

      inline FloatType compute_mean_( int ii ){
        FloatType result, tmp;
        if (centric_[ii]){
          tmp = alpha_[ii]*f_calc_[ii]/std::sqrt(2.0*epsilon_[ii]*beta_[ii]);
          result = std::exp(-tmp*tmp)*
            std::sqrt(2.0*epsilon_[ii]*beta_[ii]/scitbx::constants::pi);
          result = result + alpha_[ii]*f_calc_[ii]*scitbx::math::erf(tmp);
        }
        else {
          FloatType x;
          FloatType eb = epsilon_[ii]*beta_[ii];
          MMTBX_ASSERT(eb != 0);
          x = alpha_[ii]*f_calc_[ii]*alpha_[ii]*f_calc_[ii]/(2.0*eb);

          result = (epsilon_[ii]*beta_[ii] +
                    alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii])*
            scitbx::math::bessel::ei0(x);

          result+= alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii]*
            scitbx::math::bessel::ei1(x);

          result = result*0.5*std::sqrt(scitbx::constants::pi/(eb));
        }
        return result;

      }

      inline FloatType compute_sigma_( int ii )
      {
        // sigma^2 = epsilon_[ii] beta_[ii] + alpha_[ii]*f_calc_[ii] - mean_f_obs^2
        FloatType result;
        result = epsilon_[ii]*beta_[ii] + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii];
        result = result - mean_fo_[ii]*mean_fo_[ii];
        MMTBX_ASSERT(result >= 0);
        result = std::sqrt( result );
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
      scitbx::af::shared<FloatType> f_obs_fst_der_;
      scitbx::af::shared<FloatType> f_obs_snd_der_;

      //scitbx::af::shared<FloatType> f_obs_expected_log_likelihood_;


      scitbx::af::shared<FloatType> mean_fo_;
      scitbx::af::shared<FloatType> std_fo_;
  };




  template<typename FloatType>
  class sigmaa_estimator
  {
  public:
    sigmaa_estimator( scitbx::af::const_ref<FloatType> const& e_obs,
                      scitbx::af::const_ref<FloatType> const& e_calc,
                      scitbx::af::const_ref<bool>      const& centric,
                      scitbx::af::const_ref<FloatType> const& d_star_cubed,
                      FloatType const& width )
      :
      current_h_(-1),
      width_(width),
      eps_(1e-5),
      sigmaa_upper_bound_(1-eps_*eps_),
      log_two_(std::log(2.0)),
      log_two_minus_log_pi_( std::log(2.0)-std::log(scitbx::constants::pi) )

      {
        FloatType tmp_obs, tmp_calc;
        SCITBX_ASSERT(width > 0);
        SCITBX_ASSERT( e_obs.size() == e_calc.size() );
        SCITBX_ASSERT( e_obs.size() == centric.size() );
        SCITBX_ASSERT( e_obs.size() == d_star_cubed.size() );
        for (int ii=0;ii<e_obs.size();ii++){
          SCITBX_ASSERT( e_obs[ii]>= 0);
          SCITBX_ASSERT( e_calc[ii]>= 0);
          tmp_obs  = e_obs[ii];
          tmp_calc = e_calc[ii];
          if (tmp_obs  < 1e-3){ tmp_obs  = 1e-3; }
          if (tmp_calc < 1e-3){ tmp_calc = 1e-3; }
          shrd_e_obs_.push_back( tmp_obs );
          shrd_e_calc_.push_back( tmp_calc );
          shrd_centric_.push_back( centric[ii] ) ;
          shrd_d_star_cubed_.push_back( d_star_cubed[ii] );
          shrd_distance_.push_back(0.0);
        }

        e_obs_  = shrd_e_obs_.const_ref();
        e_calc_ = shrd_e_calc_.const_ref();
        centric_= shrd_centric_.const_ref();
        d_star_cubed_ = shrd_d_star_cubed_.const_ref();
        distance_ = shrd_distance_.const_ref();
       }

      void
      update_current_h(FloatType const& h)
      {
        FloatType delta = h-current_h_;
        if (delta<0){
          delta = -delta;
        }
        if (delta > eps_){
          current_h_ = h;
          recompute_distance();
        }
      }

      FloatType
      sum_weights(FloatType const& h)
      {
        update_current_h(h);
        return scitbx::af::sum(distance_);
      }

      FloatType target(FloatType const& h,
                       FloatType const& sigmaa
                      )
      {
        update_current_h(h);
        FloatType result=0;
        for (int ii=0;ii<e_obs_.size();ii++){
          result += distance_[ii] * compute_single_target(
            ii, std::min(sigmaa, sigmaa_upper_bound_));
        }
        return(result);
      }

      FloatType
      dtarget(FloatType const& h,
              FloatType const& sigmaa )
      {
        update_current_h(h);
        FloatType result=0;
        for (int ii=0;ii<e_obs_.size();ii++){
          result += distance_[ii] * compute_single_dtarget(
            ii, std::min(sigmaa, sigmaa_upper_bound_));
        }
        return(result);
      }

      scitbx::af::tiny<FloatType,2> target_and_gradient( FloatType const& h,
                                                         FloatType const& sigmaa )
        {
          scitbx::af::tiny<FloatType,2> result, tmp;
          result[0]=0;
          result[1]=0;
          update_current_h(h);
          for (int ii=0;ii<e_obs_.size();ii++){
            tmp = compute_single_target_and_gradient(ii,sigmaa);
            result[0]+=distance_[ii] *tmp[0];
            result[1]+=distance_[ii] *tmp[1];
          }
          return result;
        }


  protected:
      //--------------------------


      scitbx::af::tiny<FloatType,2> compute_single_target_and_gradient(int const& ii,
                                                                       FloatType sigmaa)
        {
          if (centric_[ii]){
            return( target_and_gradient_single_centric( ii, sigmaa) );
          } else {
            return( target_and_gradient_single_acentric( ii, sigmaa) );
          }
        }

      scitbx::af::tiny<FloatType,2> target_and_gradient_single_acentric(int const& ii,
                                                                        FloatType sigmaa)
        {
          scitbx::af::tiny<FloatType,2> result;
          FloatType target,grad,x,a1,a2,a3,tmp0,tmp1,tmp2,eos,ecs,eoc;
          //precompute some stuff
          eos=e_obs_[ii]*e_obs_[ii];
          ecs=e_calc_[ii]*e_calc_[ii];
          eoc=e_obs_[ii]*e_calc_[ii];
          tmp0=sigmaa*sigmaa;
          if (tmp0>=1-1e-8){
            tmp0 = 1.0-1e-8;
          }
          tmp1=1.0/(1.0-tmp0);
          tmp2=tmp1*tmp1;
          x = sigmaa*2.0*eoc*tmp1;
          // target part
          target  = log_two_ + std::log(e_obs_[ii]) + std::log( tmp1 );
          target += -(eos + tmp0*ecs)*tmp1;
          target += scitbx::math::bessel::ln_of_i0(x);
          // gradient part
          a1  = 2.0*sigmaa*tmp1;
          a2 = -2.0*sigmaa*(eos + ecs)*tmp2;
          a3 = -2.0*eoc*(1.0+tmp0)*scitbx::math::bessel::i1_over_i0(-x)*tmp2;
          grad = a1+a2+a3;
          result[0]=target;
          result[1]=grad;
          return (result);
        }

      scitbx::af::tiny<FloatType,2> target_and_gradient_single_centric(int const& ii,
                                                                       FloatType sigmaa)
        {
          scitbx::af::tiny<FloatType,2> result;
          FloatType target,grad,x,a1,a2,a3,tmp0,tmp1,tmp2,eos,ecs,eoc;
          //precompute some stuff
          eos=e_obs_[ii]*e_obs_[ii];
          ecs=e_calc_[ii]*e_calc_[ii];
          eoc=e_obs_[ii]*e_calc_[ii];
          tmp0=sigmaa*sigmaa;
          if (tmp0>=1){
            tmp0 = 0.99999;
          }
          tmp1=1.0/(1.0-tmp0);
          tmp2=tmp1*tmp1;
          x = sigmaa*eoc*tmp1;
          a1= (log_two_minus_log_pi_ + std::log(tmp1) );
          a2=-(eos+tmp0*ecs)*tmp1;
          if (x>40){
            // a simple approximation to avoid problems
            a3 = x*0.999921 - 0.65543;
          } else {
            a3 = std::log( std::cosh( x ) );
          }
          target = (a1+a2)/2.0+a3;
          // gradient part
          a1=sigmaa*tmp1;
          a2=-sigmaa*(eos+ecs)*tmp2;
          a3=eoc*(1+tmp0)*std::tanh(x)*tmp2;
          grad=a1+a2+a3;
          result[0]=target;
          result[1]=grad;
          return(result);
        }






      inline FloatType compute_single_target(int const& ii,
                                             FloatType const& sigmaa)
        {
          if (centric_[ii]){
            return( target_single_centric( ii, sigmaa) );
          } else {
            return( target_single_acentric( ii, sigmaa) );
          }
        };

      inline FloatType target_single_acentric(int const& ii,
                                              FloatType const& sigmaa
                                             )
      {
        FloatType result,x, tmp;
        tmp = 1.0-sigmaa*sigmaa;
        if (tmp<=0){
          tmp = 1e-8;
        }
        x = sigmaa*2.0*e_obs_[ii]*e_calc_[ii]/(1.0-sigmaa*sigmaa);
        result  = std::log(2.0)+std::log(e_obs_[ii])-std::log( 1.0-sigmaa*sigmaa );
        result += -(e_obs_[ii]*e_obs_[ii] + sigmaa*sigmaa*e_calc_[ii]*e_calc_[ii])/tmp;
        result += scitbx::math::bessel::ln_of_i0(x);
        return(result);
      }

      inline FloatType target_single_centric(int const& ii,
                                             FloatType const& sigmaa)
      {
        FloatType result,x,a1,a2,a3,tmp;
        tmp = 1.0-sigmaa*sigmaa;
        if (tmp<=0){
          tmp = 1e-8;
        }
        x = sigmaa*e_obs_[ii]*e_calc_[ii]/tmp;
        a1= (std::log(2.0) - std::log(scitbx::constants::pi) - std::log(1.0-sigmaa*sigmaa) )/2.0;
        a2=-(e_obs_[ii]*e_obs_[ii] +
             sigmaa*sigmaa*e_calc_[ii]*e_calc_[ii])/(2.0*tmp);
        if (x>40){
          // a simple approximation to avoid problems
          a3 = x*0.999921 - 0.65543;
        } else {
          a3 = std::log( std::cosh( x ) );
        }
        result = a1+a2+a3;
        return(result);
      }
    //-----------------------------

      FloatType compute_single_dtarget(int const& ii,
                                              FloatType const& sigmaa)
      {
        if (centric_[ii]){
          return( dtarget_single_centric(ii,sigmaa) );
        } else {
          return( dtarget_single_acentric(ii,sigmaa) );
        }
      }


      FloatType dtarget_single_acentric( int const& ii,
                                         FloatType const& sigmaa)
      {
        FloatType result=0, x,a1,a2,a3,tmp;
        tmp = 1.0-sigmaa*sigmaa;
        if (tmp<=0){
          tmp = 1e-8;
        }
        x=2.0*sigmaa*e_obs_[ii]*e_calc_[ii]/(-tmp);
        a1  = 2.0*sigmaa/tmp;
        a2 = -2.0*sigmaa*(e_obs_[ii]*e_obs_[ii] +
                          e_calc_[ii]*e_calc_[ii])/(tmp*tmp);
        a3 = -2.0*e_obs_[ii]*e_calc_[ii]*(1+sigmaa*sigmaa)*
          scitbx::math::bessel::i1_over_i0(x)/(tmp*tmp);
        result = a1+a2+a3;
        return(result);
      }


      FloatType dtarget_single_centric( int const& ii,
                                        FloatType const& sigmaa)
      {
        FloatType result, x,a1,a2,a3,tmp;
        tmp = 1.0-sigmaa*sigmaa;
        if (tmp<=0){
          tmp = 1e-8;
        }
        x = e_obs_[ii]*e_calc_[ii]*sigmaa/tmp;
        a1=sigmaa/(1-sigmaa*sigmaa);;
        a2=-sigmaa*(e_obs_[ii]*e_obs_[ii] +
                    e_calc_[ii]*e_calc_[ii])/(tmp*tmp);
        a3= e_obs_[ii]*e_calc_[ii]*(1+sigmaa*sigmaa)*std::tanh(x)/(tmp*tmp);
        result =  a1+a2+a3;
        return(result);
      }

     //-----------------------------



      inline void recompute_distance()
        {
          FloatType result;
          for (int ii=0; ii<distance_.size();ii++){
            result = current_h_- d_star_cubed_[ii];
            result = result/width_;
            result = std::exp(-result*result/2.0);
            shrd_distance_[ii]=result;
          }
        }


      scitbx::af::const_ref<FloatType> e_obs_;
      scitbx::af::const_ref<FloatType> e_calc_;
      scitbx::af::const_ref<FloatType> centric_;
      scitbx::af::const_ref<FloatType> d_star_cubed_;
      scitbx::af::const_ref<FloatType> distance_;

      scitbx::af::shared<FloatType> shrd_e_obs_;
      scitbx::af::shared<FloatType> shrd_e_calc_;
      scitbx::af::shared<FloatType> shrd_centric_;
      scitbx::af::shared<FloatType> shrd_d_star_cubed_;
      scitbx::af::shared<FloatType> shrd_distance_;

      FloatType current_h_;
      FloatType width_;
      FloatType eps_;
      FloatType sigmaa_upper_bound_;
      FloatType log_two_;
      FloatType log_two_minus_log_pi_;
  };










}}}
#endif // MMTBX_SCALING_OUTLIER_H
