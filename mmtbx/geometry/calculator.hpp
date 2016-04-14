#include <mmtbx/geometry/sphere_surface_sampling.hpp>
#include <mmtbx/geometry/indexing.hpp>
#include <mmtbx/geometry/asa.hpp>
#include <mmtbx/geometry/containment.hpp>

#include <boost/range.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_const.hpp>

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
      indexer_.add( index, coordinate_accessor_[ index ] );
    }
  }

  radius_type operator()(size_type const& index) const
  {
    namespace mxg = mmtbx::geometry;

    coordinate_type const& centre( coordinate_accessor_[ index ] );
    radius_type radius( radius_accessor_[ index ] + probe_ );

    typedef typename indexer_type::range_type close_objects_range_type;
    close_objects_range_type cor( indexer_.close_to( centre ) );

    typedef utility::Sphere< coordinate_type > sphere_type;
    typedef mxg::containment::Checker< sphere_type, mxg::containment::PurePythagorean< false > >
      checker_type;
    typedef typename boost::range_iterator< close_objects_range_type >::type
      close_objects_range_iterator;

    checker_type checker;

    for (
      close_objects_range_iterator it = boost::begin( cor );
      it != boost::end( cor );
      ++it
      )
    {
      size_type cid( *it );

      if ( cid == index )
      {
        continue;
      }

      coordinate_type const& o_centre( coordinate_accessor_[ cid ] );
      radius_type o_radius( radius_accessor_[ cid ] + probe_ );


      if ( overlap_between_spheres( centre, radius, o_centre, o_radius ) )
      {
        checker.add( sphere_type( o_centre, o_radius ) );
      }
    }

    std::ptrdiff_t dist = boost::distance(
        boost::adaptors::filter(
          boost::adaptors::transform(
            sampling_.points(),
            mxg::asa::Transform< coordinate_type >( centre, radius )
            ),
          checker
          )
        );

    return sampling_.unit_area() * radius * radius * boost::distance(
      boost::adaptors::filter(
        boost::adaptors::transform(
          sampling_.points(),
          mxg::asa::Transform< coordinate_type >( centre, radius )
          ),
        checker
        )
      );
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

} // namespace calculator

} // namespace asa

} // namespace geometry

} // namespace mmtbx
