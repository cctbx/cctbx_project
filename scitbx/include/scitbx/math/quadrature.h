#ifndef SCITBX_MATH_QUADRATURE
#define SCITBX_MATH_QUADRATURE

#include <iostream>

namespace scitbx{
namespace math{
namespace quadrature{

  // here we quickly determine the roots of a hermite polynome
  template <typename FloatType>
  class gauss_hermite_engine{
  public:
    gauss_hermite_engine(int const& n)
    {
      SCITBX_ASSERT( n < 30 ); // only stable/correct nodes are found for polynomes under 30
      p_const_ = 0.7511255444649425;
      conv_limit_ = 1.0e-13;
      n_=n;

      long m=long( (n+1)/2 ); // the number of zeros we are looking for, excluding the origin
      FloatType step,x,w;
      bool odd=true;
      if ( long(n) == std::floor(FloatType(n)/2.0)*2 ){
        odd=false;
      }

      // This last point estimate is from NR
      // it is not really needed, nor used, as shown below.
      // The roots for the hermite polynomes seem to be correct and unique
      // up to about polynhomes of 124, further has not been investigated.

      /*
      last = std::sqrt( 2.0*n_+1.0 )-1.85575*std::pow( 2.0*n_+1,-0.166667);
      last = refine( last );
      */

      if (odd){
        x_.push_back( 0.0 );
        w = f(0.0)[1];
        w_.push_back(2.0/(w*w));
        // for odd n, the first root (besides the origin) can be guessed by
        step = 2.0/std::sqrt(static_cast<double>(n_));
        step  = refine( step );
        x_.push_back( step );
        w = f(step)[1];
        w_.push_back(2.0/(w*w));
        for (int ii=2;ii<m;ii++){
          x = refine( x_[ii-1]+1.3*step );
          w = f(x)[1];
          x_.push_back( x );
          w_.push_back( 2.0/(w*w) );
          step = x_[ii] - x_[ii-1];
        }
        for (int ii=1;ii<m;ii++){
          x_.push_back( -x_[ii] );
          w_.push_back(  w_[ii] );
        }
      }

      if (!odd){
        step = 1.0/std::sqrt(static_cast<double>(n));
        step  = refine( step );
        w = f(step)[1];
        x_.push_back(step);
        w_.push_back(2.0/(w*w));
        if (n_>2){
          step = step + 2.0*step;
          step  = refine(step);
          x_.push_back(step);
          w = f(step)[1];
          w_.push_back(2.0/(w*w));

          for (int ii=2;ii<m;ii++){
            step = x_[ii-1] - x_[ii-2];
            x = refine( x_[ii-1]  + 1.3*step  );
            w = f(x)[1];
            x_.push_back( x );
            w_.push_back( 2.0/(w*w) );
          }
          for (int ii=0;ii<m;ii++){
            x_.push_back( -x_[ii] );
            w_.push_back(  w_[ii] );
          }
        }
      }

      fillit();

    }


    FloatType refine(FloatType const& z)
    {
      long max_iter = 100000, count=0;
      FloatType max_step = 1.0/( 2.0*std::sqrt(FloatType(n_)) ) ; // a mamixum step size
      FloatType delta=1, z1,z0,dd, ss;
      std::vector<FloatType> fdf;
      z1 = z;
      ss=1.0;
      while (delta > conv_limit_ ){
        z0 = z1;
        fdf = f( z1 );
        dd = fdf[0]/( fdf[1]+1e-13 );
        if ( std::fabs(dd) >= max_step ){
          if (dd<0){
            ss=-1.0;
          } else {
            ss = 1.0;
          }
          dd = max_step*ss;
        }

        z1 = z0 - dd; // newton step
        delta = std::fabs( z1-z0 );
        count++;
        if (count > max_iter ){ // this actually never happens and is just a fail safe
          delta = 0.0;
        }
      }

      return (z1);
    }


    std::vector<FloatType> f(FloatType const& z)
    {
      FloatType p1,p2,p3,pp;
      p1=p_const_;
      p2=0.0;
      for (int ii=0;ii<n_;ii++){
        p3=p2;
        p2=p1;
        p1=z*sqrt(2.0/(ii+1.0))*p2-std::sqrt(( (ii+1)-1.0)/(ii+1))*p3;
      }
      std::vector<FloatType> result;
      result.push_back( p1 ); // the function value
      // we need a derivative as well for fixed point iterations
      pp=std::sqrt(2.0*n_)*p2;
      result.push_back( pp );
      return( result );
    }

    void
    fillit()
    {
      for (int ii=0;ii<n_;ii++){
        w_exp_x_squared_.push_back( w_[ii]*std::exp( x_[ii]*x_[ii] ) );
      }
    }


    scitbx::af::shared<FloatType>
    x(){ return(x_); }

    scitbx::af::shared<FloatType>
    w(){ return(w_); }

    scitbx::af::shared<FloatType>
    w_exp_x_squared(){ return( w_exp_x_squared_ ); }

  protected:
    scitbx::af::shared<FloatType> x_;
    scitbx::af::shared<FloatType> w_;
    scitbx::af::shared<FloatType> w_exp_x_squared_;
    long n_;
    FloatType p_const_, conv_limit_;

  };



}}}
#endif
