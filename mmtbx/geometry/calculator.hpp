#ifndef MMTBX_GEOMETRY_ASA_CALCULATOR_H
#define MMTBX_GEOMETRY_ASA_CALCULATOR_H

#include <mmtbx/geometry/sphere_surface_sampling.hpp>
#include <mmtbx/geometry/indexing.hpp>
#include <mmtbx/geometry/asa.hpp>
#include <mmtbx/geometry/containment.hpp>

#include <boost/range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>

#include <stdexcept>

namespace mmtbx
{

namespace geometry
{

namespace asa
{

namespace calculator
{

namespace utility
{

template< typename Array, typename Converter >
class TransformedArray : private Converter
{
public:
  typedef Array array_type;
  typedef Converter converter_type;

  typedef typename boost::remove_reference< array_type >::type::size_type size_type;
  typedef typename converter_type::value_type return_value_type;
  typedef typename boost::remove_const<
    typename boost::remove_reference< return_value_type >::type
    >::type
    value_type;

private:
  array_type const array_;

public:
  TransformedArray(array_type const& array) : array_( array )
  {};

  return_value_type operator [](size_type const& index) const
  {
    return converter_type::operator ()( array_[ index ] );
  }

  size_type size() const
  {
    return array_.size();
  }
};

template< typename Value, typename Size = std::size_t >
class SingleValueArray
{
public:
  typedef Value value_type;
  typedef Size size_type;
  typedef value_type const& return_value_type;

private:
  value_type value_;
  size_type size_;

public:
  SingleValueArray(value_type const& value, size_type const& size)
  : value_( value ), size_( size )
  {};

  return_value_type operator [](size_type const& index) const
  {
    return value_;
  }

  size_type size() const
  {
    return size_;
  }
};

template< typename Vector >
class Sphere
{
public:
  typedef Vector vector_type;
  typedef typename vector_type::value_type value_type;

private:
  vector_type centre_;
  value_type radius_sq_;

public:
  Sphere(vector_type const& centre, value_type const& radius)
  : centre_( centre ), radius_sq_( radius * radius )
  {};

  inline const vector_type& centre() const
  {
    return centre_;
  }

  inline const value_type& radius_sq() const
  {
    return radius_sq_;
  }
};

} // namespace utility

template< typename XyzAccess, typename RadiusAccess, typename Discrete = int >
class SimpleCalculator
{
public:
  typedef XyzAccess coordinate_access_type;
  typedef RadiusAccess radius_access_type;
  typedef Discrete discrete_type;

  typedef typename boost::remove_reference< coordinate_access_type >::type::size_type size_type;
  typedef typename boost::remove_reference< coordinate_access_type >::type::value_type coordinate_type;
  typedef typename boost::remove_reference< radius_access_type >::type::value_type radius_type;

  typedef mmtbx::geometry::sphere_surface_sampling::GoldenSpiral< coordinate_type > sampling_type;
  typedef mmtbx::geometry::indexing::Hash< size_type, coordinate_type, discrete_type > indexer_type;

private:
  coordinate_access_type const coordinate_accessor_;
  radius_access_type const radius_accessor_;

  radius_type probe_;
  sampling_type sampling_;
  indexer_type indexer_;

public:
  SimpleCalculator(
    coordinate_access_type const& coords,
    radius_access_type const& radii,
    radius_type probe = 1.4,
    std::size_t sampling_point_count = 960,
    radius_type cubesize = 7.0,
    int margin = 1
    )
    : coordinate_accessor_( coords ),
      radius_accessor_( radii ),
      probe_( probe ),
      sampling_( sampling_point_count ),
      indexer_(
        typename indexer_type::voxelizer_type(
          0 < coords.size() ? coords[0] : coordinate_type( 0, 0, 0 ),
          coordinate_type( cubesize, cubesize, cubesize )
          ),
        margin
        )
  {
    size_type size( coords.size() );
    assert ( size == radii.size() );

    for (size_type index = 0; index < size; ++index )
    {
      if ( 0 < radius_accessor_[ index ] )
      {
        indexer_.add( index, coordinate_accessor_[ index ] );
      }
    }
  }

  std::ptrdiff_t accessible_points(size_type const& index) const
  {
    radius_type radius( radius_accessor_[ index ] );

    if ( radius < 0 )
    {
      throw std::runtime_error( "Requested position set to IGNORE (negative radius)" );
    }

    radius += probe_;

    namespace mxg = mmtbx::geometry;

    typedef utility::Sphere< coordinate_type > sphere_type;
    typedef mxg::containment::Checker< sphere_type, mxg::containment::PurePythagorean< false > >
      checker_type;
    checker_type checker;

    typedef typename indexer_type::range_type close_objects_range_type;
    typedef typename boost::range_iterator< close_objects_range_type >::type
      close_objects_range_iterator;
    coordinate_type centre( coordinate_accessor_[ index ] );
    close_objects_range_type cor( indexer_.close_to( centre ) );

    size_type o_index;
    coordinate_type o_centre;
    radius_type o_radius;

    for (
      close_objects_range_iterator it = boost::begin( cor );
      it != boost::end( cor );
      ++it
      )
    {
      o_index = *it;
      o_radius = radius_accessor_[ o_index ];

      if ( ( o_index == index ) || ( o_radius < 0 ) )
      {
        continue;
      }

      o_centre = coordinate_accessor_[ o_index ];
      o_radius += probe_;

      if ( overlap_between_spheres( centre, radius, o_centre, o_radius ) )
      {
        checker.add( sphere_type( o_centre, o_radius ) );
      }
    }

    return boost::distance(
      boost::adaptors::filter(
        boost::adaptors::transform(
          sampling_.points(),
          mxg::asa::Transform< coordinate_type >( centre, radius )
          ),
        checker
        )
      );
  }

  radius_type accessible_surface_area(size_type const& index) const
  {
    std::ptrdiff_t dist = accessible_points( index );
    radius_type radius( radius_accessor_[ index ] + probe_ );
    return sampling_.unit_area() * radius * radius * dist;
  }

  bool is_overlapping_sphere(coordinate_type const& centre, radius_type const& radius) const
  {
    namespace mxg = mmtbx::geometry;

    typedef typename indexer_type::range_type close_objects_range_type;
    typedef typename boost::range_iterator< close_objects_range_type >::type
      close_objects_range_iterator;
    close_objects_range_type cor( indexer_.approx_within_sphere( centre, radius ) );

    size_type o_index;
    coordinate_type o_centre;
    radius_type o_radius;

    for (
      close_objects_range_iterator it = boost::begin( cor );
      it != boost::end( cor );
      ++it
      )
    {
      o_index = *it;
      o_radius = radius_accessor_[ o_index ];

      if ( o_radius < 0 )
      {
        continue;
      }

      o_centre = coordinate_accessor_[ o_index ];

      if ( overlap_between_spheres( centre, radius, o_centre, o_radius ) )
      {
        return true;
      }
    }

    return false;
  }

private:
  bool overlap_between_spheres(
    coordinate_type left_c,
    radius_type left_r,
    coordinate_type right_c,
    radius_type right_r
  ) const
  {
    radius_type sum_radii = left_r + right_r;
    return ( left_c - right_c ).length_sq() < sum_radii * sum_radii;
  }
};

template< typename XyzAccess, typename RadiusType = double, typename Discrete = int >
class ConstRadiusCalculator
{
public:
  typedef XyzAccess coordinate_access_type;
  typedef RadiusType radius_type;
  typedef Discrete discrete_type;

  typedef typename boost::remove_reference< coordinate_access_type >::type::value_type coordinate_type;
  typedef typename boost::remove_reference< coordinate_access_type >::type::size_type size_type;
  typedef utility::SingleValueArray< radius_type, size_type > radius_access_type;
  typedef SimpleCalculator<
    coordinate_access_type,
    radius_access_type,
    discrete_type
    > calculator_type;

private:
  calculator_type calculator_;

public:
  ConstRadiusCalculator(
    const coordinate_access_type& coords,
    radius_type const& radius,
    radius_type probe = 1.4,
    std::size_t sampling_point_count = 960,
    radius_type cubesize = 7.0,
    int margin = 1
    )
    : calculator_(
        coords,
        radius_access_type( radius, coords.size() ),
        probe,
        sampling_point_count,
        cubesize,
        margin
        )
  {}

  std::ptrdiff_t accessible_points(size_type const& index) const
  {
    return calculator_.accessible_points( index );
  }

  radius_type accessible_surface_area(size_type const& index) const
  {
    return calculator_.accessible_surface_area( index );
  }

  bool is_overlapping_sphere(coordinate_type const& centre, radius_type const& radius) const
  {
    return calculator_.is_overlapping_sphere( centre, radius );
  }
};

} // namespace calculator

} // namespace asa

} // namespace geometry

} // namespace mmtbx

#endif // MMTBX_GEOMETRY_ASA_CALCULATOR_H
