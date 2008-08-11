#ifndef SCITBX_MATH_HALTON_H
#define SCITBX_MATH_HALTON_H

#include <scitbx/array_family/shared.h>
#include <scitbx/error.h>
#include <vector>

namespace scitbx{ namespace math {
namespace halton {

  /*
   * An implemenmtaionm of the halton sequence
   * for multidimensional numerical integration.
   * The sequence is known only to work for low dimensional systems,
   * so the high dimension limit is set to 6.
   * only limited python interfaces, as this routine will
   * be mostly used from C++
   */
  template <typename FloatType>
  class halton{
  public:
    halton( int const& dimension)
    {
      SCITBX_ASSERT(dimension>0);
      SCITBX_ASSERT(dimension<=6); // the sequence is not suited for
                                  // high dimensions. The limit 6 is okai ...
      FloatType tmp_primes[] = {2,3,5,7,11,13};
      for (int ii=0;ii<dimension;ii++){
        bases_.push_back( tmp_primes[ii] );
      }
      dimso_ = dimension;
    }

    FloatType nth_given_base(int const& base,
                             int const& n)
    {
      FloatType b=FloatType(base);
      FloatType i=n,half=1.0/b;
      FloatType h=0,digit;
      while (i>0){
        digit=int(i)%base;
        h += digit*half;
        i  = (i-digit)/b;
        half /= b;
      }
      return(h);
    }

    std::vector<FloatType> nth_all(int const& n)
    {
      std::vector<FloatType> result;
      for (int ii=0;ii<dimso_;ii++){
        result.push_back( nth_given_base( bases_[ii], n) );
      }
      return( result  );
    }

    scitbx::af::shared< std::vector<FloatType> >
    sequence( long const& n )
    {
      scitbx::af::shared< std::vector< FloatType > > result;
      for (long ii=0;ii<n;ii++){
        result.push_back( nth_all(ii) );
      }
      return(result);
    }


  protected:
    std::vector<FloatType> bases_;
    int dimso_;
  };


  template<typename FloatType>
  class square_halton_sampling
  {
  public:
    square_halton_sampling(
      FloatType const& low_x,
      FloatType const& high_x,
      FloatType const& low_y,
      FloatType const& high_y):
      low_x_(low_x),
      low_y_(low_y),
      high_x_(high_x),
      high_y_(high_y),
      halton_(2),
      state_(0)
      {}

    scitbx::af::tiny<FloatType,2>
    next()
    {
      std::vector<FloatType> tmp;
      tmp = halton_.nth_all(state_)  ;
      state_++;
      scitbx::af::tiny<FloatType,2> tmp_2(
        low_x_ + (high_x_ - low_x_)*tmp[0],
        low_y_ + (high_y_ - low_y_)*tmp[1]);
      return( tmp_2 );
    }
    scitbx::af::tiny<FloatType,2>
    start()
    {
      state_=0;
      return ( next() );

    }

    int state()
    {
      return(state_);
    }

    void set_state(int n)
    {
      SCITBX_ASSERT(n>=0);
      state_=n;
    }

  protected:
    int state_;
    FloatType low_x_;
    FloatType low_y_;
    FloatType high_x_;
    FloatType high_y_;
    halton<FloatType> halton_;
  };

}}} //scitbx::math::halton

#endif // SCITBX_MATH_HALTON_H
