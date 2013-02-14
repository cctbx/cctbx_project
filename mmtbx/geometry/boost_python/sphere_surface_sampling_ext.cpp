#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>

#include <scitbx/vec3.h>

#include <mmtbx/geometry/sphere_surface_sampling.hpp>
#include <boost_adaptbx/boost_range_python.hpp>


namespace mmtbx
{
namespace geometry
{
namespace sphere_surface_sampling
{
namespace
{

template < typename Vector >
struct sphere_surface_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    boost_adaptbx::python::generic_range_wrapper<
      typename GoldenSpiral< Vector >::storage_type
      >::wrap( "points_range" );

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
        make_function(
          &GoldenSpiral< Vector >::points,
          return_internal_reference<>()
          )
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

