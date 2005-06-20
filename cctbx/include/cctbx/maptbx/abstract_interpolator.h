#ifndef CCTBX_MAPTBX_ABSTRACT_INTERPOLATOR_H
#define CCTBX_MAPTBX_ABSTRACT_INTERPOLATOR_H
// Done by Jacob Smith && Erik McKee

#include <exception>
#include <cctbx/coordinates.h>
#include <cctbx/maptbx/eight_point_interpolation.h>
#if !(defined(__linux__) && defined(__GNUC__) \
 && __GNUC__ == 2 && __GNUC_MINOR__ == 96)
#include <limits>
#endif

namespace cctbx { namespace maptbx { namespace abstract {

// This file provides a mechanism for abstracting both crystallographic
// and non-crystallographic interpolation of real-space electron density
// [fourier] maps. By "crystallographic" or "non-crystallogrpaphic" we
// mean whether or not it is "symmetry aware" (==crystallographic).
//
// This is a policy based, magic-type-reinstantiation, visitor template
// class. Please see your references.
//
// interface: this policy defines the partial template specialization
//            which defines the pure abstract virtual base class
// non_crystallographic: this policy defines a partial template
//            specialization which is non-crystallographic, and therefore
//            not symmetry aware
// crystallographic<cartesian>: this policy defines a partial template
//            specialization which is crystallographic (symmetry aware),
//            specifically the map has been "fractionalized" and is thus
//            a compact cartesian subset
// crystallographic<fractional>: this policy defines a partial template
//            specialization which is crystallographic (symmetry aware),
//            but the map is in "fractional" space and thus the cartesian
//            input coordinates must be converted to fractional coordinates
//            before being used
//
// The user should use a interpolator<FloatType[,void]> ONLY, DO NOT USE
// THE OTHER TYPES!!! the interpolator<FloatType[,void]> will do all the
// stuff you need (copy construction, assignment, construction, etc. etc.).
// The interpolator<FloatType[,void]> will decide what sort of back-end
// object to use based upon the constructor/available data.
//
// If you want to extend this, here is what you must do:
// define a policy:
//     struct my_policy {};
//
// create a partial template specialization that inherits from interpolator<interface,FloatType>:
//     template < typename FloatType >
//     class interpolator<my_policy,FloatType> : public interpolate<interface,FloatType> {
//       // override virtual functions
//     };
//
// create a parametrized overload of the get_*_interpolator functor:
//     template < typename FloatType >
//     struct get_my_policy_interpolator {
//       static
//       interpolator<FloatType>
//       interpolator ( ..., my_policy const& mp = my_policy() ) {
//         return interpolator<FloatType>( interpolator<my_policy,FloatType>(...) );
//       }
//     };
//

struct there_is_no_interpolator_backend : public std::exception {};

#if BOOST_WORKAROUND(__EDG_VERSION__, == 238)
class edg_workaround_cartesian { friend class cartesian; };
class edg_workaround_fractional { friend class fractional; };
#endif

struct interface {};
struct cartesian {};
struct fractional {};
template < typename T=void > struct crystallographic {};
template < typename T=void > struct non_crystallographic {};

template < typename T, typename U=void > class interpolator;

template < typename FloatType > class interpolator<interface,FloatType> {
public:
  virtual ~ interpolator () {}
  virtual FloatType interpolate ( scitbx::vec3<FloatType> const& ) const = 0;
  virtual interpolator<interface,FloatType>* clone () const = 0;
};

template < typename FloatType >
class interpolator<non_crystallographic<>,FloatType>
  : public interpolator<interface,FloatType> {
public:
  interpolator ( af::versa<FloatType, af::flex_grid<> > const& Map,
                 scitbx::mat3<FloatType> const& gmtx,
                 bool aoob=false,
                 FloatType const& oobsv=0 )
  : map_(Map)
  , gridding_matrix_(gmtx)
  , allow_out_of_bounds_(aoob)
  , out_of_bounds_substitute_value_(oobsv) {
    this->map_c_ref_ = this->map_.const_ref();
  }
  virtual ~ interpolator () {}
  virtual interpolator<interface,FloatType>* clone () const {
    return new interpolator<non_crystallographic<>,FloatType>(this->map_,
      this->gridding_matrix_,
      this->allow_out_of_bounds_,
      this->out_of_bounds_substitute_value_);
  }
  virtual FloatType interpolate ( scitbx::vec3<FloatType> const& site ) const {
    return cctbx::maptbx::non_crystallographic_eight_point_interpolation(
      this->map_c_ref_,
      this->gridding_matrix_,
      site,
      this->allow_out_of_bounds_,
      this->out_of_bounds_substitute_value_);
  }
private:
  af::versa<FloatType, af::flex_grid<> > map_;
  af::const_ref<FloatType, af::flex_grid<> > map_c_ref_;
  scitbx::mat3<FloatType> gridding_matrix_;
  bool allow_out_of_bounds_;
  FloatType out_of_bounds_substitute_value_;
};

template < typename FloatType >
class interpolator<crystallographic<fractional>,FloatType>
  : public interpolator<interface,FloatType> {
public:
  interpolator ( af::versa<FloatType, af::flex_grid<> > const& Map )
  : map_(Map)
  , map_c_ref_(
      af::const_ref<double, af::c_grid_padded<3> >(
        Map.begin(),
        af::c_grid_padded<3>(Map.accessor()))) {
  }
  virtual ~ interpolator () {}
  virtual interpolator<interface,FloatType>* clone () const {
    return new interpolator<crystallographic<fractional>,FloatType>(this->map_);
  }
  virtual FloatType interpolate ( scitbx::vec3<FloatType> const& site ) const {
    // theoretically this works ... perhaps Ralf has some input?
    // the idea is that fractional <==> cartesian
    cctbx::fractional<FloatType> frac(site);
    return cctbx::maptbx::eight_point_interpolation(this->map_c_ref_,frac);
  }
private:
  af::versa<FloatType, af::flex_grid<> > map_;
  af::const_ref<FloatType, af::c_grid_padded<3> > map_c_ref_;
};

template < typename FloatType >
class interpolator<crystallographic<cartesian>,FloatType>
  : public interpolator<interface,FloatType> {
public:
  interpolator ( af::versa<FloatType, af::flex_grid<> > const& Map,
    scitbx::mat3<FloatType> const& fmtx )
  : map_(Map)
  , map_c_ref_(
      af::const_ref<double, af::c_grid_padded<3> >(
        Map.begin(),
        af::c_grid_padded<3>(Map.accessor())))
  , fractionalization_matrix_(fmtx) {
  }
  virtual ~ interpolator () {}
  virtual interpolator<interface,FloatType>* clone () const {
    return new interpolator<crystallographic<cartesian>,FloatType>(this->map_,this->fractionalization_matrix_);
  }
  virtual FloatType interpolate ( scitbx::vec3<FloatType> const& site ) const {
    cctbx::fractional<FloatType> frac = this->fractionalization_matrix_ * site;
    return cctbx::maptbx::eight_point_interpolation(this->map_c_ref_,frac);
  }
private:
  af::versa<FloatType, af::flex_grid<> > map_;
  af::const_ref<FloatType, af::c_grid_padded<3> > map_c_ref_;
  scitbx::mat3<FloatType> fractionalization_matrix_;
};

template < typename FloatType >
class interpolator<FloatType,void> {
public:
  // STL container compatibility
  interpolator () : interpolator_(0) {}
  // copy c-tor for copying-goodness
  interpolator ( interpolator<FloatType> const& I ) {
    this->interpolator_ = 0;
    this->copy(I);  // <-- this is a clever bit
  }
  // assignment goodness
  interpolator<FloatType>& operator = ( interpolator<FloatType> const& I ) {
    this->copy(I);
    return *this;
  }
  // so that you can derive from the class
  virtual ~ interpolator () {
    this->clear();
  }
  // normal constructor of goodness
  template < typename T >
  interpolator ( interpolator<T,FloatType> const& I ) {
    this->interpolator_ = I.clone();
  }
  // the interface itself!!!
  FloatType interpolate ( scitbx::vec3<FloatType> const& site ) const {
    if ( this->interpolator_ != 0 )
      return this->interpolator_->interpolate(site);
    throw there_is_no_interpolator_backend();
  }
  af::shared<FloatType> interpolate ( af::const_ref<scitbx::vec3<FloatType> > const& sites ) const {
    af::shared<FloatType> result(sites.size());
    std::size_t sz = result.size();
    for ( std::size_t i=0; i<sz; ++i )
      result[i] = this->interpolate( sites[i] );
    return result;
  }
  FloatType operator [] ( scitbx::vec3<FloatType> const& site ) const {
    return this->interpolate( site );
  }
  af::shared<FloatType> operator [] ( af::const_ref<scitbx::vec3<FloatType> > const& sites ) const {
    return this->interpolate( sites );
  }
  bool can_interpolate () const { return this->interpolator_!=0; }
protected:
  void copy ( interpolator<FloatType> const& I ) {
    this->clear();
    if ( I.interpolator_ != 0 )
      this->interpolator_ = I.interpolator_->clone();
  }
  void clear () {
    if ( this->interpolator_ != 0 )
      delete this->interpolator_;
    this->interpolator_ = 0;
  }
private:
  interpolator<interface,FloatType> *interpolator_;
};

template < typename FloatType >
struct get_non_crystallographic {

  static
  interpolator<FloatType>
  interpolator ( af::versa<FloatType, af::flex_grid<> > const& Map,
                 scitbx::mat3<FloatType> const& gmtx,
                 bool aoob=false,
                 FloatType const& oobsv=0 ) {
    return abstract::interpolator<FloatType>(
              abstract::interpolator<non_crystallographic<>,FloatType>(Map,gmtx,aoob,oobsv) );
  }

};

template < typename FloatType >
struct get_fractional_crystallographic {

  static
  interpolator<FloatType>
  interpolator ( af::versa<FloatType, af::flex_grid<> > const& Map ) {
    return abstract::interpolator<FloatType>( abstract::interpolator<crystallographic<fractional>,FloatType>(Map) );
  }

};

template < typename FloatType >
struct get_cartesian_crystallographic {

  static
  interpolator<FloatType>
  interpolator ( af::versa<FloatType, af::flex_grid<> > const& Map,
                     scitbx::mat3<FloatType> const& fmtx ) {
    return abstract::interpolator<FloatType>(
            abstract::interpolator<crystallographic<cartesian>,FloatType>(Map,fmtx) );
  }

};

}}}// end cctbx::maptbx::abstract

#endif//CCTBX_MAPTBX_ABSTRACT_INTERPOLATOR_H
