#include <vector>

#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>
#include <boost/python/stl_iterator.hpp>

#include <scitbx/vec3.h>

#include <mmtbx/geometry/containment.hpp>
#include <mmtbx/geometry/overlap.hpp>
#include <mmtbx/geometry/primitive.hpp>

namespace mmtbx
{

namespace geometry
{

namespace containment
{

namespace
{

template < typename Vector >
struct containment_wrappers
{
  static void wrap()
  {
    using namespace boost::python;

    object containment_module(
      handle<>( borrowed( PyImport_AddModule( "volume.containment" ) ) )
      );
    scope().attr( "containment" ) = containment_module;
    scope containment_scope = containment_module;

    class_< PointInOriginSphere >( "point_in_origin_sphere", no_init )
      .def( init<>() )
      .def(
        "__call__",
        &PointInOriginSphere::operator ()< Vector >,
        ( arg( "point" ), arg( "radius_sq" ) )
        )
      ;

    class_< PointInOriginDiamond >( "point_in_origin_diamond", no_init )
      .def( init<>() )
      .def(
        "__call__",
        &PointInOriginDiamond::operator ()< Vector >,
        ( arg( "point" ), arg( "radius" ) )
        )
      ;

    class_< PointInOriginCube >( "point_in_origin_cube", no_init )
      .def( init<>() )
      .def(
        "__call__",
        &PointInOriginCube::operator ()< Vector >,
        ( arg( "point" ), arg( "radius" ) )
        )
      ;
  }
};

void init_module()
{
  containment_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace containment

namespace overlap
{

namespace
{

template < typename Vector >
struct overlap_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    object overlap_module(
      handle<>( borrowed( PyImport_AddModule( "volume.overlap" ) ) )
      );
    scope().attr( "overlap" ) = overlap_module;
    scope overlap_scope = overlap_module;

    using namespace mmtbx::geometry::primitive;

    class_< BetweenSpheres >( "between_spheres", no_init )
      .def( init<>() )
      .def(
        "__call__",
        &BetweenSpheres::operator ()< Sphere< Vector > >,
        ( arg( "left" ), arg( "right" ) )
        )
      ;
    class_< BetweenBoxes >( "between_boxes", no_init )
      .def( init<>() )
      .def(
        "__call__",
        &BetweenBoxes::operator ()< Box< Vector > >,
        ( arg( "left" ), arg( "right" ) )
        )
      ;
  }
};

void init_module()
{
  overlap_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace overlap
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_volume_ext)
{
  mmtbx::geometry::containment::init_module();
  mmtbx::geometry::overlap::init_module();
}

