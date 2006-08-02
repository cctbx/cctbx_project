#ifndef SCITBX_MATH_HALTON_H
#define SCITBX_MATH_HALTON_H

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
      SCITBX_ASSERT(dimension<6); // the sequence is not suited for
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



}}} //scitbx::math::halton


#endif // SCITBX_MATH_HALTON_H
