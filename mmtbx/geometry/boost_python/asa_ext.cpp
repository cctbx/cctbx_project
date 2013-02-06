#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

#include <scitbx/vec3.h>
#include <scitbx/math/boost_python/cartesian_product_fixed_size.hpp>

#include <mmtbx/geometry/asa.hpp>
#include <mmtbx/geometry/indexing.hpp>
#include <mmtbx/geometry/overlap.hpp>
#include <mmtbx/geometry/containment.hpp>
#include <mmtbx/geometry/sphere_surface_sampling.hpp>

#include <string>

namespace mmtbx
{
namespace geometry
{
namespace asa
{
namespace
{

template < typename Vector >
struct asa_wrappers
{
  static void wrap()
  {
    using namespace boost::python;

    typedef typename primitive::Traits< Vector >::vector_type vector_type;
    typedef typename primitive::Traits< Vector >::value_type value_type;

    class_< Sphere< Vector >, bases< primitive::Sphere< Vector > > >(
      "sphere",
      no_init
      )
      .def(
        init< const vector_type&, const value_type&, size_t >(
          ( arg( "centre" ), arg( "radius" ), arg( "index" ) )
          )
        )
      .add_property(
        "index",
        make_function(
          &Sphere< Vector >::index,
          return_value_policy< copy_const_reference >()
          )
        )
      .def( self == self )
      .def( self != self )
      .def(
        "create",
        Sphere< Vector >::create,
        ( arg( "centre" ), arg( "radius" ) )
        )
      .staticmethod( "create" )
      ;
  }
};

void init_module()
{
  asa_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace asa

namespace indexing
{

namespace
{

template < typename Vector >
struct indexing_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    object indexing_module(
      handle<>( borrowed( PyImport_AddModule( "asa.indexing" ) ) )
      );
    scope().attr( "indexing" ) = indexing_module;
    scope indexing_scope = indexing_module;

    typedef asa::Sphere< Vector > sphere_type;
    typedef Linear< sphere_type, overlap::BetweenSpheres > linear_spheres_type;

    scitbx::math::cartesian_product::python::iterated_range_wrappers<
      typename linear_spheres_type::const_iterator
      >::wrap( "linear_spheres_overlapping_objects_range" );

    class_< linear_spheres_type >( "linear_spheres", no_init )
      .def( init<>() )
      .def( "add", &linear_spheres_type::add, arg( "object" ) )
      .def(
        "overlapping_with",
        &linear_spheres_type::overlapping_with,
        arg( "object" )
        )
      .def( "__len__", &linear_spheres_type::size )
      ;
    
    typedef Voxelizer< Vector, scitbx::vec3< int > > voxelizer_type;

    class_< voxelizer_type >( "voxelizer", no_init )
      .def(
        init< const Vector&, const Vector& >( ( arg( "base" ), arg( "step" ) ) )
        )
      .def( "__call__", &voxelizer_type::operator (), arg( "vector" ) )
      ;
  }
};

void init_module()
{
  indexing_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace indexing

namespace containment
{

namespace
{

template< typename Checker >
void
add_neighbours_from_list(Checker& checker, boost::python::object neighbours)
{
  typedef typename Checker::neighbour_type neighbour_type;
  checker.add(
    boost::python::stl_input_iterator< neighbour_type >( neighbours ),
    boost::python::stl_input_iterator< neighbour_type >()
    );
}

template< typename Checker, typename Range >
void
add_neighbours_from_range(Checker& checker, const Range& neighbours)
{
  checker.add( neighbours.begin, neighbours.end );
}

template< typename Checker >
typename Checker::storage_type::const_iterator
checker_neighbours_begin(const Checker& checker)
{
  return checker.neighbours().begin();
}

template< typename Checker >
typename Checker::storage_type::const_iterator
checker_neighbours_end(const Checker& checker)
{
  return checker.neighbours().end();
}

template < typename Vector >
struct accessibility_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    object accessibility_module(
      handle<>( borrowed( PyImport_AddModule( "asa.accessibility" ) ) )
      );
    scope().attr( "accessibility" ) = accessibility_module;
    scope accessibility_scope = accessibility_module;

    typedef asa::Sphere< Vector > sphere_type;

    create_wrappings< sphere_type, PurePythagorean< false > >( "pythagorean" );
    create_wrappings< sphere_type, DiamondPrefilter< false > >(
      "diamond_prefilter"
      );
  }

  template< typename Sphere, typename Algorithm >
  static void create_wrappings(const std::string& name)
  {
    typedef Checker< Sphere, Algorithm > checker_type;
    typedef typename checker_type::functor_type functor_type;
    typedef typename Sphere::vector_type vector_type;

    using namespace sphere_surface_sampling;
    typedef typename GoldenSpiral< vector_type >::range_type
      transformed_range_type;
    typedef FilterResult< transformed_range_type, functor_type >
      golden_spiral_filter_helper;
    using namespace scitbx::math::cartesian_product::python;
    iterated_range_wrappers<
      typename golden_spiral_filter_helper::filter_iterator
      >::wrap( ( name + "_filtered_transformed_range" ).c_str() );

    typedef indexing::FilterHelper<
      typename std::vector< Sphere >::const_iterator,
      typename indexing::OverlapEqualityFilter<
        Sphere,
        overlap::BetweenSpheres
        >
      >
      filter_helper;

    using namespace boost::python;

    class_< checker_type >( ( name + "_checker" ).c_str(), no_init )
      .def( init<>() )
      .def(
        "add_from_list",
        add_neighbours_from_list< checker_type >,
        arg( "neighbours" )
        )
      .def(
        "add_from_range",
        add_neighbours_from_range<
          checker_type,
          typename filter_helper::filter_range_type >,
        arg( "neighbours" )
        )
      .def(
        "neighbours", 
        range(
          checker_neighbours_begin< checker_type >,
          checker_neighbours_end< checker_type >
          )
        )
      .def(
        "filter",
        &checker_type::template filter< transformed_range_type >,
        with_custodian_and_ward_postcall< 0, 2 >(),
        arg( "range" )
        )
      .def(
        "is_selected",
        &checker_type::operator (),
        arg( "point" )
        )
      ;
  }
};

void init_module()
{
  accessibility_wrappers< scitbx::vec3< double > >::wrap();
}

} // namespace <anonymous>
} // namespace containment
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_asa_ext)
{
  mmtbx::geometry::asa::init_module();
  mmtbx::geometry::indexing::init_module();
  mmtbx::geometry::containment::init_module();
}

