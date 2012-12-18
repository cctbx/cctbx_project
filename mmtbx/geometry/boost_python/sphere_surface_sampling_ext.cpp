#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/iterator.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/vec3.h>

#include <mmtbx/geometry/sphere_surface_sampling.hpp>
#include <scitbx/math/boost_python/cartesian_product_fixed_size.hpp>

namespace mmtbx
{
namespace geometry
{
namespace sphere_surface_sampling
{
namespace
{

template< typename Vector >
typename GoldenSpiral< Vector >::const_storage_iterator
sampling_points_begin(const GoldenSpiral< Vector >& sampling)
{
  return sampling.points().begin();
}

template< typename Vector >
typename GoldenSpiral< Vector >::const_storage_iterator
sampling_points_end(const GoldenSpiral< Vector >& sampling)
{
  return sampling.points().end();
}

template < typename Vector >
struct sphere_surface_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    scitbx::math::cartesian_product::python::iterated_range_wrappers<
      typename GoldenSpiral< Vector >::const_iterator
      >::wrap( "point_range" );

    class_< GoldenSpiral< Vector > >( "golden_spiral", no_init )
      .def(
        init< size_t >( ( arg( "count" ) ) )
          )
      .add_property(
        "unit_area",
        make_function(
          &GoldenSpiral< Vector >::unit_area,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "count",
        make_function(
          &GoldenSpiral< Vector >::count,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "points",
        range(
          sampling_points_begin< Vector >,
          sampling_points_end< Vector >
          )
        )
      .def(
        "transformed",
        &GoldenSpiral< Vector >::transformed,
        with_custodian_and_ward_postcall< 0, 1 >(),
        ( arg( "centre" ), arg( "radius" ) )
        )
      ;
  }
};

void init_module()
{
  sphere_surface_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace sphere_surface_sampling
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_sphere_surface_sampling_ext)
{
  mmtbx::geometry::sphere_surface_sampling::init_module();
}

