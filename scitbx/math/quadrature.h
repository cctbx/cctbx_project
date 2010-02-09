#ifndef SCITBX_MATH_QUADRATURE_H
#define SCITBX_MATH_QUADRATURE_H

#include <scitbx/array_family/tiny.h>
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/ref.h>
#include <vector>

namespace scitbx{ namespace math{ namespace quadrature{

   /***********************************************************************
   *     2D integration, in the plane, with Gaussian weight
   *     - nine_twentyone_1012
   *     - five_nine_1110
   *     - five_nine_1001
   *     - seven_twelve_0120
   *
   *

   *  This file deal with approximation the following integral:
   *  Q = \int_{-\infty}^{\infty} exp(-x^2) f(x)
   *  The approximation is carried out using a quadrature (1D)
   *  or cubature (2D).
   *  The integral is effectively approximated by
   *  Q \approx \sum_{i=1}^{N} w_i f(x_i)
   *  the elements w_i are the weights and x_i are the coordinates
   *  on which the function has to be evaluated.
   *  The weights and coordinates are givfen by the implemented classes
   *
   *  Note the following rule when changing variables:
   *  \int exp(-x^2)f(x) dx = \int J*exp(-(z-z0)/(2s^2))f(z) dz
   *  with
   *  J = dx/dz; ( in this case, J = Sqrt[2] s)
   *
   */

  // This class generates cubature point sets
  // see http://www.cs.kuleuven.be/~nines/research/ecf/
  //
  // the names of the classes are as follows:
  // order_numberofpoints_rule
  // This is in lnie with the tables given on the above website.
  // only a selected number (4) of cubatures have been implemented.
  // using 21, 12,9 and 9 points.
  // initial testing showed that the 9 point cubature did not show a large accuracy
  // improvement over using a standard multiplicative Gauss Hermite quadrature
  //
  // the functions coord() and weight() are for python tests only.
  //

  template <typename FloatType>
  class nine_twentyone_1012{
    public:
      nine_twentyone_1012()
      {
        scitbx::af::tiny<FloatType,2> tmp_1( 0,0 );
        FloatType w1=1.0471975511965977;
        x_.push_back(tmp_1);
        w_.push_back(w1);

        scitbx::af::tiny<FloatType,2> tmp_2(1.5381890013208515,1.5381890013208515);
        FloatType w2=0.012238016879947463;
        expand(tmp_2,w2,false,false);

        scitbx::af::tiny<FloatType,2> tmp_3(0.43091398228826897,1.0403183802565386);
        FloatType w3=0.24426215416421333;
        expand(tmp_3,w3,false,true);

        scitbx::af::tiny<FloatType,2> tmp_4(0.54093734825794375,2.1069972930282899);
        FloatType w4=0.011418225194962370;
        expand(tmp_4,w4,false,true);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > coords()
      {
        return( x_ );
      }

      scitbx::af::tiny<FloatType,2> coord(int const& ii)
      {
        SCITBX_ASSERT( ii<21 );
        return( x_[ii] );
      }


      scitbx::af::shared< FloatType > weights()
      {
        return( w_ );
      }

      FloatType weight( int const& ii )
      {
        SCITBX_ASSERT( ii<21 );
        return( w_[ii] );
      }

    protected:

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > plus_minus( scitbx::af::tiny<FloatType,2> const& gen,
                                                                      bool const& last_is_zero )
      {
        scitbx::af::shared<  scitbx::af::tiny<FloatType,2> > result;

        scitbx::af::tiny<FloatType,2> s1(gen[0], gen[1]);
        scitbx::af::tiny<FloatType,2> s2(-gen[0], gen[1]);
        result.push_back( s1 );
        result.push_back( s2 );

        if (!last_is_zero){
          scitbx::af::tiny<FloatType,2> s3(gen[0], -gen[1]);
          scitbx::af::tiny<FloatType,2> s4(-gen[0], -gen[1]);
          result.push_back( s3 );
          result.push_back( s4 );
        }

        return(result);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > swap_expand(
         scitbx::af::shared< scitbx::af::tiny<FloatType,2> > const& gen)
      {
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > result;
        int size = gen.size();
        for (int ii=0;ii<size;ii++){
          scitbx::af::tiny<FloatType,2> s1(gen[ii][1],gen[ii][0]);
          result.push_back( s1 );
          result.push_back( gen[ii] );
        }
        return( result );
      }

      void expand( scitbx::af::tiny<FloatType,2> const& gen,
                   FloatType const& w,
                   bool const& last_one_is_zero,
                   bool const& swapit )
      {
        // first, get the sign thing going
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > tmp;
        tmp = plus_minus( gen, last_one_is_zero );
        if (swapit){
          // now do the swap expand
          tmp = swap_expand( tmp );
        }
        // now push back the whole lot
        for (int ii=0;ii<tmp.size();ii++){
          x_.push_back( tmp[ii] );
          w_.push_back( w );
        }

      }


      scitbx::af::shared< scitbx::af::tiny<FloatType,2> >  x_;
      scitbx::af::shared< FloatType > w_;


  };




  template <typename FloatType>
  class five_nine_1110{
    public:
      five_nine_1110()
      {
        scitbx::af::tiny<FloatType,2> tmp_1(0,0);
        FloatType w1=1.5707963267948966;
        x_.push_back( tmp_1 );
        w_.push_back( w1 );

        scitbx::af::tiny<FloatType,2> tmp_2(1,1);
        FloatType w2=0.19634954084936207;
        expand( tmp_2,w2,false,false);

        scitbx::af::tiny<FloatType,2> tmp_3(1.4142135623730950, 0 );
        FloatType w3=0.19634954084936207;
        expand( tmp_3,w3,true,true);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > coords()
      {
        return( x_ );
      }

      scitbx::af::tiny<FloatType,2> coord(int const& ii)
      {
        SCITBX_ASSERT( ii <= 8 );
        return( x_[ii] );
      }


      scitbx::af::shared< FloatType > weights()
      {
        return( w_ );
      }

      FloatType weight( int const& ii )
      {
        SCITBX_ASSERT( ii <= 8 );
        return( w_[ii] );
      }

    protected:

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > plus_minus( scitbx::af::tiny<FloatType,2> const& gen,
                                                                      bool const& last_is_zero )
      {
        scitbx::af::shared<  scitbx::af::tiny<FloatType,2> > result;

        scitbx::af::tiny<FloatType,2> s1(gen[0], gen[1]);
        scitbx::af::tiny<FloatType,2> s2(-gen[0], gen[1]);
        result.push_back( s1 );
        result.push_back( s2 );

        if (!last_is_zero){
          scitbx::af::tiny<FloatType,2> s3(gen[0], -gen[1]);
          scitbx::af::tiny<FloatType,2> s4(-gen[0], -gen[1]);
          result.push_back( s3 );
          result.push_back( s4 );
        }

        return(result);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > swap_expand(
         scitbx::af::shared< scitbx::af::tiny<FloatType,2> > const& gen)
      {
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > result;
        int size = gen.size();
        for (int ii=0;ii<size;ii++){
          scitbx::af::tiny<FloatType,2> s1(gen[ii][1],gen[ii][0]);
          result.push_back( s1 );
          result.push_back( gen[ii] );
        }
        return( result );
      }

      void expand( scitbx::af::tiny<FloatType,2> const& gen,
                   FloatType const& w,
                   bool const& last_one_is_zero,
                   bool const& swapit )
      {
        // first, get the sign thing going
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > tmp;
        tmp = plus_minus( gen, last_one_is_zero );
        if (swapit){
          // now do the swap expand
          tmp = swap_expand( tmp );
        }
        // now push back the whole lot
        for (int ii=0;ii<tmp.size();ii++){
          x_.push_back( tmp[ii] );
          w_.push_back( w );
        }

      }


      scitbx::af::shared< scitbx::af::tiny<FloatType,2> >  x_;
      scitbx::af::shared< FloatType > w_;


  };



  template <typename FloatType>
  class five_nine_1001{
    public:
      five_nine_1001()
      {
        // These values have been taken from the mentioned webpage
        // a ful reference is :  A.H. Stroud, Approximate calculation
        // of multiple integrals, Prentice-Hall, Englewood Cliffs, N.J., 1971.
        scitbx::af::tiny<FloatType,2> tmp_1(0,0);
        FloatType w1=1.5707963267948966;
        x_.push_back( tmp_1 );
        w_.push_back( w1 );

        scitbx::af::tiny<FloatType,2> tmp_2(1.3065629648763765,0.54119610014619698 );
        FloatType w2=0.19634954084936207;
        expand( tmp_2,w2,false,true);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > coords()
      {
        return( x_ );
      }

      scitbx::af::tiny<FloatType,2> coord(int const& ii)
      {
        SCITBX_ASSERT( ii <= 8 );
        return( x_[ii] );
      }


      scitbx::af::shared< FloatType > weights()
      {
        return( w_ );
      }

      FloatType weight( int const& ii )
      {
        SCITBX_ASSERT( ii <= 8 );
        return( w_[ii] );
      }

    protected:

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > plus_minus( scitbx::af::tiny<FloatType,2> const& gen,
                                                                      bool const& last_is_zero )
      {
        scitbx::af::shared<  scitbx::af::tiny<FloatType,2> > result;

        scitbx::af::tiny<FloatType,2> s1(gen[0], gen[1]);
        scitbx::af::tiny<FloatType,2> s2(-gen[0], gen[1]);
        result.push_back( s1 );
        result.push_back( s2 );

        if (!last_is_zero){
          scitbx::af::tiny<FloatType,2> s3(gen[0], -gen[1]);
          scitbx::af::tiny<FloatType,2> s4(-gen[0], -gen[1]);
          result.push_back( s3 );
          result.push_back( s4 );
        }

        return(result);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > swap_expand(
         scitbx::af::shared< scitbx::af::tiny<FloatType,2> > const& gen)
      {
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > result;
        int size = gen.size();
        for (int ii=0;ii<size;ii++){
          scitbx::af::tiny<FloatType,2> s1(gen[ii][1],gen[ii][0]);
          result.push_back( s1 );
          result.push_back( gen[ii] );
        }
        return( result );
      }

      void expand( scitbx::af::tiny<FloatType,2> const& gen,
                   FloatType const& w,
                   bool const& last_one_is_zero,
                   bool const& swapit )
      {
        // first, get the sign thing going
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > tmp;
        tmp = plus_minus( gen, last_one_is_zero );
        if (swapit){
          // now do the swap expand
          tmp = swap_expand( tmp );
        }
        // now push back the whole lot
        for (int ii=0;ii<tmp.size();ii++){
          x_.push_back( tmp[ii] );
          w_.push_back( w );
        }

      }


      scitbx::af::shared< scitbx::af::tiny<FloatType,2> >  x_;
      scitbx::af::shared< FloatType > w_;


  };




  template <typename FloatType>
  class seven_twelve_0120{
    public:
      seven_twelve_0120()
      {
        // These values have been taken from the mentioned webpage
        // a ful reference is :  A.H. Stroud, Approximate calculation
        // of multiple integrals, Prentice-Hall, Englewood Cliffs, N.J., 1971.
        scitbx::af::tiny<FloatType,2> tmp_1(1.7320508075688772, 0);
        FloatType w1=0.087266462599716478;
        expand( tmp_1,w1,true,true);

        scitbx::af::tiny<FloatType,2> tmp_2(0.53523313465963489, 0.53523313465963489);
        FloatType w2=0.66127983844512042;
        expand( tmp_2,w2,false,false);

        scitbx::af::tiny<FloatType,2> tmp_3(1.4012585384440735, 1.4012585384440735);
        FloatType w3=0.036851862352611409;
        expand( tmp_3,w3,false,false);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > coords()
      {
        return( x_ );
      }

      scitbx::af::tiny<FloatType,2> coord(int const& ii)
      {
        SCITBX_ASSERT( ii <= 11 );
        return( x_[ii] );
      }


      scitbx::af::shared< FloatType > weights()
      {
        return( w_ );
      }

      FloatType weight( int const& ii )
      {
        SCITBX_ASSERT( ii <= 11 );
        return( w_[ii] );
      }

    protected:

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > plus_minus( scitbx::af::tiny<FloatType,2> const& gen,
                                                                      bool const& last_is_zero )
      {
        scitbx::af::shared<  scitbx::af::tiny<FloatType,2> > result;

        scitbx::af::tiny<FloatType,2> s1(gen[0], gen[1]);
        scitbx::af::tiny<FloatType,2> s2(-gen[0], gen[1]);
        result.push_back( s1 );
        result.push_back( s2 );

        if (!last_is_zero){
          scitbx::af::tiny<FloatType,2> s3(gen[0], -gen[1]);
          scitbx::af::tiny<FloatType,2> s4(-gen[0], -gen[1]);
          result.push_back( s3 );
          result.push_back( s4 );
        }

        return(result);
      }

      scitbx::af::shared< scitbx::af::tiny<FloatType,2> > swap_expand(
         scitbx::af::shared< scitbx::af::tiny<FloatType,2> > const& gen)
      {
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > result;
        int size = gen.size();
        for (int ii=0;ii<size;ii++){
          scitbx::af::tiny<FloatType,2> s1(gen[ii][1],gen[ii][0]);
          result.push_back( s1 );
          result.push_back( gen[ii] );
        }
        return( result );
      }

      void expand( scitbx::af::tiny<FloatType,2> const& gen,
                   FloatType const& w,
                   bool const& last_one_is_zero,
                   bool const& swapit )
      {
        // first, get the sign thing going
        scitbx::af::shared< scitbx::af::tiny<FloatType,2> > tmp;
        tmp = plus_minus( gen, last_one_is_zero );
        if (swapit){
          // now do the swap expand
          tmp = swap_expand( tmp );
        }
        // now push back the whole lot
        for (int ii=0;ii<tmp.size();ii++){
          x_.push_back( tmp[ii] );
          w_.push_back( w );
        }

      }


      scitbx::af::shared< scitbx::af::tiny<FloatType,2> >  x_;
      scitbx::af::shared< FloatType > w_;


  };



  /**********************************************************
   *    INtegration over a sphere
   *    implemented:
   *     - des.3.17.4 :   17 point spherical 3-design
   *     - des.3.24.7 :   24 point spherical 7-design
   *     - des.3.40.8 :   40 point spherical 8-design
   *     - des.3.60.10:   60 point spherical 10-design
   *     - des.3.240.21: 240 point spherical 21-design
   *
   *********************************************************/







  /********************************************************
   *
   *          1 Dimensional integration
   *    - Gauss legendre [-1,1] unit weights
   *    - Gauss hermite [-infty, infty], exp(-x*x) weights
   *
   ********************************************************/

  // here we determine roots of legendre polynomes needed for 1d
  // integration over the interval [-1,1]
  template<typename FloatType>
  class gauss_legendre_engine{
  public:
    gauss_legendre_engine(int const& n)
    {
       SCITBX_ASSERT( n < 96 );
       SCITBX_ASSERT( n > 1 );
       n_ = n;
       conv_limit_ = 1e-13;
       max_count_ = 1000; // 1000 iterations should really be enough

       FloatType start_guess=1.0-1e-5;
       FloatType x;
       for (int ii=0;ii<int( (n+1)/2 ) ;ii++){
         x = refine( start_guess );
           x_.push_back( x );
           w_.push_back( f(x)[2] );
         if ( std::fabs(x) > conv_limit_ ){
           x_.push_back( -x );
           w_.push_back( f(x)[2] );
         }
       }
    }

    // using the Newton algorithm with implicit deflating
    // a.k.a. The Durant Kerner Formula. See:
    // http://digilander.libero.it/foxes/poly/Zeros_orthogonal_polynomials.htm
    //
    FloatType refine( FloatType const& z )
    {
      FloatType delta=100,zn,zold,inflator,ratio;
      std::vector<FloatType> fdf;
      zn=z;
      int count=0;
      while ( delta > conv_limit_ ){
        zold = zn;
        inflator=0;
        for (int ii=0;ii<x_.size();ii++){
          inflator += 1.0/(zn-x_[ii]);
        }
        fdf = f(zn);
        ratio = fdf[0]/(fdf[1]-fdf[0]*inflator);
        zn = zn-ratio;
        delta = std::fabs( zn - zold );
        count++;
        if (count>= max_count_){
          delta = 0;
        }
      }

      return( zn );
    }


    std::vector<FloatType> f(FloatType const& z)
    {
      FloatType p1, p2, p3;
      p1=1.0;
      p2=0.0;
      p3=0.0;
      for (int jj=0;jj<n_;jj++){
        p3=p2;
        p2=p1;
        p1=((2.0*(jj+1)-1.0)*z*p2-((jj+1)-1.0)*p3)/(jj+1.0);
      }
      std::vector<FloatType> result;
      result.push_back( p1 );
      result.push_back( n_*(z*p1-p2)/(z*z-1.0) );
      result.push_back( 2.0/( (1.0-z*z)*result[1]*result[1] ) );
      return (result);
    }

    scitbx::af::shared<FloatType> x()
    {
      return( x_ );
    }

    scitbx::af::shared<FloatType> w()
    {
      return( w_ );
    }



  protected:
    int n_;
    int max_count_;
    FloatType conv_limit_;
    scitbx::af::shared<FloatType> x_;
    scitbx::af::shared<FloatType> w_;

  };




  // here we quickly determine the roots of a hermite polynome needed for 1d integration.
  template <typename FloatType>
  class gauss_hermite_engine{
  public:
    gauss_hermite_engine(int const& n)
    {
      SCITBX_ASSERT(n > 1 );
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
        }
        for (int ii=0;ii<m;ii++){
          x_.push_back( -x_[ii] );
          w_.push_back(  w_[ii] );
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
      SCITBX_ASSERT(x_.size() == n_);
      SCITBX_ASSERT(w_.size() == n_);
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

}}} // namespace scitbx::math::quadrature

#endif // SCITBX_MATH_QUADRATURE_H
