// copyright (c) Jacob N. Smith & Erik McKee; leave this here; use at your whim
#ifndef CCTBX_MAPTBX_COORDINATE_TRANSFORMERS_HPP
#define CCTBX_MAPTBX_COORDINATE_TRANSFORMERS_HPP

#include<scitbx/mat3.h>
#include<scitbx/math/utils.h>
#include<cctbx/coordinates.h>

namespace cctbx {

namespace maptbx {

static const std::size_t dimension_3 = 3;

// converts a grid into a fractional
template < typename FloatType, typename IntType >
fractional<FloatType>
grid_fractionalization ( grid_point<IntType> const& coordinate,
  af::tiny<IntType,dimension_3> const& extents ) {
  scitbx::vec3<FloatType> result;
  for ( std::size_t i=0; i<dimension_3; ++i )
    result[i] = FloatType(coordinate[i]) / FloatType(extents[i]);
  return result;
}

// converts a fractional into a grid
template < typename IntType, typename FloatType >
scitbx::vec3<FloatType>
strange_fractional_gridization ( fractional<FloatType> const& coordinate,
    af::tiny<IntType,dimension_3> const& extents ) {
  scitbx::vec3<FloatType> result;
  for ( std::size_t i=0; i<dimension_3; ++i )
    result[i] = coordinate[i] * extents[i];
  return result;
}

template < typename IntType, typename FloatType >
grid_point<IntType>
floor_fractional_gridization ( fractional<FloatType> const& coordinate,
    af::tiny<IntType,dimension_3> const& extents ) {
  typedef scitbx::math::float_int_conversions<FloatType,IntType> fic;
  scitbx::vec3<FloatType> intermediate = strange_fractional_gridization(coordinate,extents);
  grid_point<IntType> result;
  for ( std::size_t i=0; i<dimension_3; ++i )
    result[i] = IntType( fic::ifloor( intermediate[i] ) );
  return result;
}


template < typename IntType, typename FloatType >
grid_point<IntType>
round_fractional_gridization ( fractional<FloatType> const& coordinate,
    af::tiny<IntType,dimension_3> const& extents ) {
  typedef scitbx::math::float_int_conversions<FloatType,IntType> fic;
  scitbx::vec3<FloatType> intermediate = strange_fractional_gridization(coordinate,extents);
  grid_point<IntType> result;
  for ( std::size_t i=0; i<dimension_3; ++i )
    result[i] = IntType( fic::iround( intermediate[i] ) );
  return result;
}

template <  typename FromType,
      typename ToType> struct transform;

/* mips is t3h suxx0r */
// unblobs the transform type to get the from and to types
//template < typename T > struct transformer_types;
//template < template<typename,typename> class T, typename U1, typename U2 >
//struct transformer_types< T<U1,U2> > {
//  typedef U1 from_type;
//  typedef U2 to_type;
//};
template < typename T > struct transformer_types {
  typedef typename T::from_type from_type;
  typedef typename T::to_type to_type;
};

// T->T transformations require no work
template < typename T >
struct transform<T,T> {
  typedef T from_type;
  typedef T to_type;
  transform () {}
  T operator () ( T const& t ) const {
    return t;
  }
  maptbx::transform<T,T> inverse () const {
    return maptbx::transform<T,T>();
  }
};

// expand from fractional to grid
template < typename FloatType, typename IntType >
struct transform<fractional<FloatType>,grid_point<IntType> > {
  typedef fractional<FloatType> from_type;
  typedef grid_point<IntType> to_type;
  transform () {}
  transform ( af::tiny<IntType,dimension_3> const& to_grid ) {
    this->to_grid_ = to_grid;
  }
  grid_point<IntType> operator () ( fractional<FloatType> const& coord ) const {
    return round_fractional_gridization(coord,this->to_grid_);
  }
  maptbx::transform<grid_point<IntType>,fractional<FloatType> > inverse () const {
    return maptbx::transform<grid_point<IntType>,fractional<FloatType> >(this->to_grid_);
  }
  grid_point<IntType> floor_transform ( fractional<FloatType> const& coord ) const {
    return floor_fractional_gridization(coord,this->to_grid_);
  }
  scitbx::vec3<FloatType> strange_transform ( fractional<FloatType> const& coord ) const {
    return strange_fractional_gridization(coord,this->to_grid_);
  }
  af::tiny<IntType,dimension_3> to_grid_;
};

// fractional to cartesian
template < typename FloatType >
struct transform<fractional<FloatType>,cartesian<FloatType> > {
  typedef fractional<FloatType> from_type;
  typedef cartesian<FloatType> to_type;
  transform () {}
  transform ( scitbx::mat3<FloatType> const& to_cart ) {
    this->to_cart_ = to_cart;
  }
  cartesian<FloatType> operator () ( fractional<FloatType> const& coord ) const {
    return this->to_cart_ * coord;
  }
  maptbx::transform<cartesian<FloatType>,fractional<FloatType> > inverse () const {
    return maptbx::transform<cartesian<FloatType>,fractional<FloatType> >(this->to_cart_.inverse());
  }
  scitbx::mat3<FloatType> to_cart_;
};

// cartesian to grid
template < typename FloatType, typename IntType >
struct transform<cartesian<FloatType>,grid_point<IntType> > {
  typedef cartesian<FloatType> from_type;
  typedef grid_point<IntType> to_type;

  transform () {}
  transform (
    maptbx::transform<cartesian<FloatType>,fractional<FloatType> > const& c2f,
    maptbx::transform<fractional<FloatType>,grid_point<IntType> > const& f2g ) {
    this->to_frac_ = c2f.to_frac_;
    this->to_grid_ = f2g.to_grid_;
  }
  transform ( scitbx::mat3<FloatType> const& to_frac,
    af::tiny<IntType,dimension_3> const& to_grid ) {
    this->to_frac_ = to_frac;
    this->to_grid_ = to_grid;
  }
  grid_point<IntType> operator () ( cartesian<FloatType> const& coord ) const {
    scitbx::vec3<FloatType> F = this->to_frac_*coord;
    fractional<FloatType> frac(F[0],F[1],F[2]);
    return grid_point<IntType>( round_fractional_gridization( frac,
      this->to_grid_) );
  }
  maptbx::transform<grid_point<IntType>,cartesian<FloatType> > inverse () const {
    return maptbx::transform<grid_point<IntType>,cartesian<FloatType> >(
      this->to_grid_, this->to_frac_.inverse() );
  }
  scitbx::mat3<FloatType> to_frac_;
  af::tiny<IntType,dimension_3> to_grid_;
};

// cartesian to fractional
template < typename FloatType >
struct transform<cartesian<FloatType>,fractional<FloatType> > {
  typedef cartesian<FloatType> from_type;
  typedef fractional<FloatType> to_type;

  transform () {}
  transform ( scitbx::mat3<FloatType> const& to_frac ) {
    this->to_frac_ = to_frac;
  }
  fractional<FloatType> operator () ( cartesian<FloatType> const& coord ) const {
    return this->to_frac_ * coord;
  }
  maptbx::transform<fractional<FloatType>,cartesian<FloatType> > inverse () const {
    return maptbx::transform<fractional<FloatType>,cartesian<FloatType> >(this->to_frac_.inverse());
  }
  scitbx::mat3<FloatType> to_frac_;
};

// grid to cartesian
template < typename FloatType, typename IntType >
struct transform<grid_point<IntType>,cartesian<FloatType> > {
  typedef grid_point<IntType> from_type;
  typedef cartesian<FloatType> to_type;

  transform () {}
  transform (
    maptbx::transform<grid_point<IntType>,fractional<FloatType> > const& g2f,
    maptbx::transform<fractional<FloatType>,cartesian<FloatType> > const& f2c ) {
    this->to_cart_ = f2c.to_cart_;
    this->to_frac_ = g2f.to_frac_;
  }
  transform ( af::tiny<IntType,dimension_3> const& to_frac,
    scitbx::mat3<FloatType> const& to_cart ) {
    this->to_frac_ = to_frac;
    this->to_cart_ = to_cart;
  }
  cartesian<FloatType> operator () ( grid_point<IntType> const& coord ) const {
    return cartesian<FloatType>( this->to_cart_ *
        grid_fractionalization<FloatType>(coord,this->to_frac_) );
  }
  maptbx::transform<cartesian<FloatType>,grid_point<IntType> > inverse () const {
    return maptbx::transform<cartesian<FloatType>,grid_point<IntType> >(
      this->to_cart_.inverse(), this->to_frac_ );
  }
  scitbx::mat3<FloatType> to_cart_;
  af::tiny<IntType,dimension_3> to_frac_;
};

// grid to fractional
template < typename FloatType, typename IntType >
struct transform<grid_point<IntType>,fractional<FloatType> > {
  typedef grid_point<IntType> from_type;
  typedef fractional<FloatType> to_type;

  transform () {}
  transform ( af::tiny<IntType,dimension_3> const& to_frac ) {
    this->to_frac_ = to_frac;
  }
  fractional<FloatType> operator () ( grid_point<IntType> const& coord ) const {
    return grid_fractionalization<FloatType>(coord,this->to_frac_);
  }
  maptbx::transform<fractional<FloatType>,grid_point<IntType> > inverse () const {
    return maptbx::transform<fractional<FloatType>,grid_point<IntType> >(this->to_frac_);
  }
  af::tiny<IntType,dimension_3> to_frac_;
};

}// maptbx

}// cctbx

#endif//CCTBX_MAPTBX_COORDINATE_TRANSFORMERS_HPP
