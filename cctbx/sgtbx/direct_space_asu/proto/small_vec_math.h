#ifndef SCITBX_AF_VEC3_MATH_H
#define SCITBX_AF_VEC3_MATH_H

#include <cmath>
#include <limits>

#include <boost/static_assert.hpp>
#include <boost/rational.hpp>
#include <scitbx/vec3.h>

namespace scitbx {

  namespace {
    BOOST_STATIC_ASSERT( -3/2 == -1 ); // -3/2 == -2 in python
  }

  template< typename IntType>
  inline IntType ceil(const boost::rational<IntType> &r)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    SCITBX_ASSERT( r.denominator() > 0  );
    if( r.denominator()==1 || r.numerator()==0 )
      return r.numerator();
    SCITBX_ASSERT( r.numerator() != 0 );
    IntType i = r.numerator()/r.denominator();
    if( r.numerator()>0 ) // -3/2 == -1 in C and -2 in python
      ++i;
    return i;
  }

  template< typename IntType >
  inline IntType floor(const boost::rational<IntType> &r)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    SCITBX_ASSERT( r.denominator() > 0  );
    if( r.denominator()==1 )
      return r.numerator();
    SCITBX_ASSERT( r.numerator()!=0  );
    IntType i = r.numerator()/r.denominator();
    if( r.numerator()<0 ) // -3/2 == -1 in C and -2 in python
      --i;
    return i;
  }

  template< typename IntType >
  inline IntType iround(const boost::rational<IntType> &r)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    SCITBX_ASSERT( r.denominator() > 0  );
    if( r.denominator()==1 )
      return r.numerator();
    double i = boost::rational_cast<double,IntType>(r);
    i = std::floor( i + 0.5 );
    return static_cast<IntType> ( i );
  }


  template<typename IntType >
  inline bool operator<= (const double lhs, const boost::rational<IntType> &rhs)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    SCITBX_ASSERT( rhs.denominator()>0 );
    return lhs*rhs.denominator() <= rhs.numerator();
  }

  template<typename IntType>
  inline bool operator>= (const double lhs, const boost::rational<IntType> &rhs)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    SCITBX_ASSERT( rhs.denominator()>0 );
    return lhs*rhs.denominator() >= rhs.numerator();
  }

  template<typename ReturnType, typename InputType>
    inline scitbx::vec3<ReturnType> vec3_cast(const scitbx::vec3<InputType> &i)
  {
    return scitbx::vec3<ReturnType>( static_cast<ReturnType>( i[0] ),
        static_cast<ReturnType>( i[1] ), static_cast<ReturnType>( i[2] ) );
  }

  template< typename FloatType >
  inline scitbx::vec3<FloatType> floor(const scitbx::vec3<FloatType> &v)
  {
    BOOST_STATIC_ASSERT( !std::numeric_limits< FloatType >::is_integer );
    return scitbx::vec3<FloatType>( std::floor(v[0]), std::floor(v[1]), std::floor(v[2]) );
  }

  template< typename FloatType >
  inline scitbx::vec3<FloatType> round(const scitbx::vec3<FloatType> &v)
  {
    BOOST_STATIC_ASSERT( !std::numeric_limits< FloatType >::is_integer );
    return scitbx::vec3<FloatType>( std::floor(v[0]+0.5), std::floor(v[1]+0.5), std::floor(v[2]+0.5) );
  }


  template< typename IntType >
  inline scitbx::vec3<IntType> ceil(const scitbx::vec3< boost::rational<IntType> > &rv)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    return scitbx::vec3<IntType>( scitbx::ceil(rv[0]), scitbx::ceil(rv[1]), scitbx::ceil(rv[2]) );
  }

  template< typename IntType>
  inline scitbx::vec3<IntType> floor(const scitbx::vec3< boost::rational<IntType> > &r)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    return scitbx::vec3<IntType>( scitbx::floor(r[0]), scitbx::floor(r[1]), scitbx::floor(r[2]) );
  }

  template< typename IntType>
  inline scitbx::vec3<IntType> iround(const scitbx::vec3< boost::rational<IntType> > &r)
  {
    BOOST_STATIC_ASSERT( std::numeric_limits< IntType >::is_integer );
    return scitbx::vec3<IntType>( scitbx::iround(r[0]), scitbx::iround(r[1]), scitbx::iround(r[2]) );
  }

  template<typename LhsNumType, typename RhsNumType>
  inline bool le_all(const scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    for(signed char i=0; i<3; ++i)
      if( !(lhs[i]<=rhs[i]) )
        return false;
    return true;
  }

  template<typename LhsNumType, typename RhsNumType>
  inline bool lt_all(const scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    for(signed char i=0; i<3; ++i)
      if( !(lhs[i]<rhs[i]) )
        return false;
    return true;
  }


  template< typename LhsNumType, typename RhsNumType>
  inline bool ge_all(const scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    for(signed char i=0; i<3; ++i)
      if( !(lhs[i]>=rhs[i]) )
        return false;
    return true;
  }

  template< typename LhsNumType, typename RhsNumType>
  inline bool gt_all(const scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    for(signed char i=0; i<3; ++i)
      if( !(lhs[i]>rhs[i]) )
        return false;
    return true;
  }


  template<typename LhsNumType, typename RhsNumType>
  inline bool eq_all(const scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    for(signed char i=0; i<3; ++i)
      if( lhs[i]!=rhs[i] )
        return false;
    return true;
  }


  template<typename LhsNumType, typename RhsNumType>
  inline scitbx::vec3<LhsNumType> & operator+= (scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    lhs[0] += rhs[0]; lhs[1] += rhs[1]; lhs[2] += rhs[2];
    return lhs;
  }

  template<typename LhsNumType, typename RhsNumType>
  inline scitbx::vec3<LhsNumType> operator+ (const scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    return scitbx::vec3<LhsNumType>(lhs[0] + rhs[0], lhs[1] + rhs[1], lhs[2] + rhs[2]);
  }

  template<typename LhsNumType, typename RhsNumType>
  inline scitbx::vec3<LhsNumType> & operator-= (scitbx::vec3<LhsNumType> &lhs, const scitbx::vec3<RhsNumType> &rhs)
  {
    lhs[0] -= rhs[0]; lhs[1] -= rhs[1]; lhs[2] -= rhs[2];
    return lhs;
  }


  template<typename LhsType, typename RhsType>
  inline void mul(scitbx::vec3<LhsType> &lhs, const scitbx::vec3<RhsType> &rhs)
  {
    lhs[0] *= rhs[0]; lhs[1]*=rhs[1]; lhs[2]*=rhs[2];
  }


  template<typename NumType>
  inline std::ostream & out (std::ostream &os, const scitbx::vec3<NumType> &v)
  {
    os << "(" << v[0] << ", " << v[1] << ", " << v[2] << ")";
    return os;
  }

  template<typename NumType>
  inline std::ostream & operator<< (std::ostream &os, const scitbx::vec3<NumType> &v)
  {
    return scitbx::out(os, v);
  }


  typedef boost::rational<int> rational_int;
  typedef scitbx::vec3< signed char > tiny3;
  typedef scitbx::vec3< int > int3;
  typedef scitbx::vec3< long > long3;
  typedef scitbx::vec3< double > double3;
  typedef scitbx::vec3< rational_int > rational_int3;

}

#endif
