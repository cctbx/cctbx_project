#ifndef CCTBX_MAPTBX_MAP_OUT_OF_BOUNDS_H
#define CCTBX_MAPTBX_MAP_OUT_OF_BOUNDS_H

namespace cctbx {

namespace maptbx {

struct raise;
struct substitute;
struct clamp;
struct interpolate;

// out of bounds interface to deal with out-of-bounds accesses
template <  typename DealWithItType,
      typename FloatType,
      typename IntType >
  struct out_of_bounds;

// generalized interface
template < typename FloatType, typename IntType >
struct out_of_bounds<void,FloatType,IntType> {
  typedef out_of_bounds<void,FloatType,IntType> h_oob;
  typedef chiltbx::handle::handle<h_oob> handle_type;
  virtual handle_type as_handle () const = 0;
  // everything & the kitchen sink =)
  virtual FloatType retry (
    generic_grid<void,FloatType,IntType> const&,
    transform<grid_point<IntType>, fractional<FloatType> > const&,
      fractional<FloatType> const& ) const = 0;
};

template < typename FloatType, typename IntType >
struct out_of_bounds<raise,FloatType,IntType>
: public out_of_bounds<void,FloatType,IntType> {
  typedef out_of_bounds<void,FloatType,IntType> h_oob;
  typedef chiltbx::handle::handle<h_oob> handle_type;
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  // everything & the kitchen sink =)
  virtual FloatType retry (
    generic_grid<void,FloatType,IntType> const&,
    transform<grid_point<IntType>, fractional<FloatType> > const&,
      fractional<FloatType> const& ) const {
    throw error("basic_map<T>: the coordinate is out of bounds.");
  }
};

template < typename FloatType, typename IntType >
struct out_of_bounds<clamp,FloatType,IntType>
: public out_of_bounds<void,FloatType,IntType> {
  out_of_bounds ( FloatType const& value ) : value_(value) {}
  typedef out_of_bounds<void,FloatType,IntType> h_oob;
  typedef chiltbx::handle::handle<h_oob> handle_type;
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  // everything & the kitchen sink =)
  virtual FloatType retry (
    generic_grid<void,FloatType,IntType> const&,
    transform<grid_point<IntType>, fractional<FloatType> > const&,
      fractional<FloatType> const& ) const {
    throw this->value_;
  }
  FloatType value_;
};

template < typename FloatType, typename IntType >
struct out_of_bounds<substitute,FloatType,IntType>
: public out_of_bounds<void,FloatType,IntType> {
  out_of_bounds ( FloatType const& value ) : value_(value) {}
  typedef out_of_bounds<void,FloatType,IntType> h_oob;
  typedef chiltbx::handle::handle<h_oob> handle_type;
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  // everything & the kitchen sink =)
  virtual FloatType retry (
    generic_grid<void,FloatType,IntType> const&,
    transform<grid_point<IntType>, fractional<FloatType> > const&,
      fractional<FloatType> const& ) const {
    return this->value_;
  }
  FloatType value_;
};

template < typename FloatType, typename IntType >
struct out_of_bounds<interpolate,FloatType,IntType>
: public out_of_bounds<void,FloatType,IntType> {
  typedef out_of_bounds<void,FloatType,IntType> h_oob;
  typedef chiltbx::handle::handle<h_oob> handle_type;
  virtual handle_type as_handle () const {
    return handle_type(*this);
  }
  // everything & the kitchen sink =)
  virtual FloatType retry (
    generic_grid<void,FloatType,IntType> const&,
    transform<grid_point<IntType>, fractional<FloatType> > const&,
      fractional<FloatType> const& ) const {
    throw error("basic_mapper<T>: non-trivial interpolation is not implemented.");
  }
};

}

}

#endif//CCTBX_MAPTBX_MAP_OUT_OF_BOUNDS_H
