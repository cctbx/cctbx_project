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


namespace mmtbx { namespace scaling { namespace outlier{

  // Some model based outlier criteria are computed here
  template <typename FloatType>
  class likelihood_ratio_outlier_test{
  public:
    likelihood_ratio_outlier_test(
          scitbx::af::const_ref<FloatType> const& f_obs,
          scitbx::af::const_ref<FloatType> const& s_obs,
          scitbx::af::const_ref<FloatType> const& f_calc,
          scitbx::af::const_ref<FloatType> const& epsilon,
          scitbx::af::const_ref<bool> const& centric,
          scitbx::af::const_ref<FloatType> const& alpha,
          scitbx::af::const_ref<FloatType> const& beta ):
      log_ei0_(40000) // the exponentiated besseli0 function (  exp[-x]Io[x]  )
      {
        FloatType min_beta=10.0;
        FloatType min_alpha=1e-5;//
        /*
        scitbx::math::quadrature::gauss_legendre_engine<FloatType>
          gauss_legendre_(15);
        x_ = gauss_legendre_.x();
        w_ = gauss_legendre_.w();
        */
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
          // get the stuff please
          newton_search(ii, 0.50,  1e-5, 10000);
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
        bool converged=false, use_mean=false;
        int count=0;
        start_f_bar = mean_fo_[ii]*start_fraction;
        f_bar = mean_fo_[ii]*start_fraction;
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
          if (step>0.5*start_f_bar){
            step=0.5*start_f_bar;
          }
          if (step < -0.5*start_f_bar){
            step=-0.5*start_f_bar;
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
            /*
            std::cout << " PRELIMINAIRY CONVERGENCE   "
                      << delta << " "
                      << centric_[ii] << " " << f_calc_[ii] << " "
                      << f_bar << " " << " " <<  f_bar_m1 << " " <<  f_bar_m2 << " "
                      << mean_fo_[ii] <<  " "
                      << alpha_[ii]   << " "
                      << beta_[ii] << " " << std::endl;
            */
            use_mean=true;
            converged=true;
          }
          count++;
          //std::cout << count <<" " << f_bar << " " << funct << " " << grad << "  " << -funct/grad << std::endl;

        }
        if (use_mean){
          //std::cout <<"using f_mean rather then mode"<< std::endl;
          f_bar = mean_fo_[ii];
        }
        grad  = snd_der(f_bar, ii);
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
        FloatType eb=epsilon_[ii]*beta_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x = alpha_[ii]*fo*f_calc_[ii]/(eb);
        FloatType m = std::tanh(x)*alpha_[ii]*f_calc_[ii]/(eb);;
        result = -fo/eb + m;
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
                 - (2.0/eb)
                 + (f_calc_[ii]*4.0*alpha_[ii]*alpha_[ii]/(eb*eb))*
          (1.0-m/(2.0*x)-m*m);
        return (result);
      }

      inline FloatType calc_snd_der_centric( FloatType fo, int ii )
      {
        FloatType result=0;
        FloatType eb=epsilon_[ii]*beta_[ii];
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
        FloatType eb=epsilon_[ii]*beta_[ii];
        if (fo<=1e-13){
          fo=1e-13;
        }
        FloatType x=alpha_[ii]*fo*f_calc_[ii]/eb;
        FloatType exparg = (fo*fo + alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii])/(2.0*eb);
        result = 0.5*std::log(2.0)-0.5*std::log(scitbx::constants::pi) - 0.5*std::log(eb)
          -exparg + std::log(std::cosh(x));
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
          x = alpha_[ii]*f_calc_[ii]*alpha_[ii]*f_calc_[ii]/(
            2.0*epsilon_[ii]*beta_[ii]);

          result = (epsilon_[ii]*beta_[ii] +
                    alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii])*
            scitbx::math::bessel::ei0(x);

          result+= alpha_[ii]*alpha_[ii]*f_calc_[ii]*f_calc_[ii]*
            scitbx::math::bessel::ei1(x);

          result = result*0.5*std::sqrt(
            scitbx::constants::pi/(epsilon_[ii]*beta_[ii]));

        }
        return result;

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


      mmtbx::scaling::twinning::quick_log_ei0<FloatType> log_ei0_;
      scitbx::af::shared<FloatType> mean_fo_;

  };




  template<typename FloatType>
  class sigmaa_estimator
  {
  public:
    sigmaa_estimator( scitbx::af::const_ref<FloatType> const& e_obs,     //1
                      scitbx::af::const_ref<FloatType> const& e_calc,    //2
                      scitbx::af::const_ref<bool>      const& centric,   //4
                      scitbx::af::const_ref<FloatType> const& d_star_sq, //6
                      FloatType const& width )                           //7
      :
      current_h_(0.0),
      width_(width)
      {
        SCITBX_ASSERT( e_obs.size() == e_calc.size() );
        SCITBX_ASSERT( e_obs.size() == centric.size() );
        SCITBX_ASSERT( e_obs.size() == d_star_sq.size() );
        eps_=1e-5;
        for (int ii=0;ii<e_obs.size();ii++){
          SCITBX_ASSERT( e_obs[ii]> 0);
          SCITBX_ASSERT( e_calc[ii]> 0);

          e_obs_.push_back( e_obs[ii] );
          e_calc_.push_back( e_calc[ii] );
          centric_.push_back( centric[ii] ) ;
          d_star_sq_.push_back( d_star_sq[ii] );
          distance_.push_back(0.0);

        }
      }

      FloatType target(FloatType const& h,
                       FloatType const& sigmaa
                      )
      {
        FloatType delta = h-current_h_;
        if (delta<0){
          delta = -delta;
        }
        if (delta > eps_){
          current_h_ = h;
          recompute_distance();
        }
        FloatType result=0, tmp_sigmaa=sigmaa;
        if (sigmaa>1.0-eps_*eps_){
          tmp_sigmaa = 1.0-eps_*eps_;
        }
        for (int ii=0;ii<e_obs_.size();ii++){
          result += distance_[ii]*compute_single_target(ii,tmp_sigmaa);
        }
        return(result);
      }

      FloatType
      dtarget(FloatType const& h,
              FloatType const& sigmaa )
      {
        FloatType delta = h-current_h_;
        if (delta<0){
          delta = -delta;
        }
        if (delta > eps_){
          current_h_ = h;
          recompute_distance();
        }

        FloatType result=0, tmp_sigmaa=sigmaa;
        if (sigmaa>1.0-eps_*eps_){
          tmp_sigmaa = 1.0-eps_*eps_;
        }
        for (int ii=0;ii<e_obs_.size();ii++){
          result += distance_[ii]*compute_single_dtarget(ii,tmp_sigmaa);
        }

        return(result);
      }

  protected:
      //--------------------------

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
        FloatType result,x;
        x = sigmaa*2.0*e_obs_[ii]*e_calc_[ii]/(1.0-sigmaa*sigmaa);
        result  = std::log(2.0)+std::log(e_obs_[ii])-std::log( 1.0-sigmaa*sigmaa );
        result += -(e_obs_[ii]*e_obs_[ii] + sigmaa*sigmaa*e_calc_[ii]*e_calc_[ii])/
          (1.0-sigmaa*sigmaa);
        result += scitbx::math::bessel::ln_of_i0(x);
        return(result);
      }

      inline FloatType target_single_centric(int const& ii,
                                             FloatType const& sigmaa)
      {
        FloatType result,x,a1,a2,a3;
        x = sigmaa*e_obs_[ii]*e_calc_[ii]/(1.0-sigmaa*sigmaa);
        a1= (std::log(2.0) - std::log(scitbx::constants::pi) - std::log(1.0-sigmaa*sigmaa) )/2.0;
        a2=-(e_obs_[ii]*e_obs_[ii] +
             sigmaa*sigmaa*e_calc_[ii]*e_calc_[ii])/(2.0*(1.0-sigmaa*sigmaa));
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
        FloatType result=0, x,a1,a2,a3;
        x=2.0*sigmaa*e_obs_[ii]*e_calc_[ii]/(sigmaa*sigmaa-1.0);
        a1  = 2.0*sigmaa/(1-sigmaa*sigmaa);
        a2 = -2.0*sigmaa*(e_obs_[ii]*e_obs_[ii] +
                          e_calc_[ii]*e_calc_[ii])/((1.0-sigmaa*sigmaa )*(1.0-sigmaa*sigmaa ));
        a3 = -2.0*e_obs_[ii]*e_calc_[ii]*(1+sigmaa*sigmaa)*
          scitbx::math::bessel::i1_over_i0(x)/((1.0-sigmaa*sigmaa)*(1.0-sigmaa*sigmaa));
        result = a1+a2+a3;
        return(result);
      }


      FloatType dtarget_single_centric( int const& ii,
                                        FloatType const& sigmaa)
      {
        FloatType result, x,a1,a2,a3;
        x = e_obs_[ii]*e_calc_[ii]*sigmaa/(1.0-sigmaa*sigmaa);
        a1=sigmaa/(1-sigmaa*sigmaa);;
        a2=-sigmaa*(e_obs_[ii]*e_obs_[ii] +
                    e_calc_[ii]*e_calc_[ii])/((1.0-sigmaa*sigmaa )*(1.0-sigmaa*sigmaa ));
        a3= e_obs_[ii]*e_calc_[ii]*(1+sigmaa*sigmaa)*std::tanh(x)/((1.0-sigmaa*sigmaa)*(1.0-sigmaa*sigmaa));
        result =  a1+a2+a3;
        return(result);
      }

     //-----------------------------



      inline void recompute_distance()
        {
          FloatType result;
          for (int ii=0; ii<distance_.size();ii++){
            result = current_h_- d_star_sq_[ii];
            result = result/width_;
            result = std::exp(-result*result/2.0);
            distance_[ii]=result;
          }
        }


      scitbx::af::shared<FloatType> e_obs_;
      scitbx::af::shared<FloatType> e_calc_;
      scitbx::af::shared<FloatType> centric_;
      scitbx::af::shared<FloatType> d_star_sq_;
      scitbx::af::shared<FloatType> distance_;

      FloatType current_h_;
      FloatType width_;
      FloatType eps_;
  };










}}}
#endif // MMTBX_SCALING_OUTLIER_H
