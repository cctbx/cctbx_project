#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <scitbx/vec3.h>

#include <mmtbx/geometry/primitive.hpp>

namespace mmtbx
{

namespace geometry
{

namespace primitive
{

namespace
{

template < typename Vector >
struct primitive_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    using namespace mmtbx::geometry::primitive;

    typedef typename Traits< Vector >::vector_type vector_type;
    typedef typename Traits< Vector >::value_type value_type;

    class_< Sphere< Vector > >( "sphere", no_init )
      .def(
        init< const vector_type&, const value_type& >(
          ( arg( "centre" ), arg( "radius" ) )
          )
        )
      .add_property(
        "centre",
        make_function(
          &Sphere< Vector >::centre,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "radius",
        make_function(
          &Sphere< Vector >::radius,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "radius_sq",
        make_function(
          &Sphere< Vector >::radius_sq,
          return_value_policy< copy_const_reference >()
          )
        )
      ;

    class_< Box< Vector > >( "box", no_init )
      .add_property(
        "low",
        make_function(
          &Box< Vector >::low,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "high",
        make_function(
          &Box< Vector >::high,
          return_value_policy< copy_const_reference >()
          )
        )
      .def(
        "from_corners",
        &Box< Vector >::from_corners,
        ( arg( "corner1" ), arg( "corner2" ) )
        )
      .staticmethod( "from_corners" )
      .def(
        "around_sphere",
        &Box< Vector >::around_sphere,
        ( arg( "centre" ), arg( "radius" ) )
        )
      .staticmethod( "around_sphere" )
      ;

    class_< BSphere< Vector >, bases< Sphere< Vector >, Box< Vector > > >(
      "bsphere",
      no_init
      )
      .def(
        init< const vector_type&, const value_type& >(
          ( arg( "centre" ), arg( "radius" ) )
          )
        )
      ;
  }
};

void init_module()
{
  primitive_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace primitive
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_primitive_ext)
{
  mmtbx::geometry::primitive::init_module();
}

