#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/stl_iterator.hpp>

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
      ;

    class_< IndexFilter< Vector > >( "index_filter", no_init )
      .def( init< const size_t& >( arg( "index" ) ) )
      .def( "__call__", &IndexFilter< Vector >::operator (), arg( "object" ) )
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

    typedef PrefilterHelper<
      typename linear_spheres_type::storage_type::const_iterator,
      asa::IndexFilter< Vector >,
      typename linear_spheres_type::overlap_filter_type
      >
      prefilter_helper;

    scitbx::math::cartesian_product::python::iterated_range_wrappers<
      typename prefilter_helper::filter_range_type::iterator_type
      >::wrap( "linear_spheres_overlapping_index_prefiltered_objects_range" );

    class_< linear_spheres_type >( "linear_spheres", no_init )
      .def( init<>() )
      .def( "add", &linear_spheres_type::add, arg( "object" ) )
      .def(
        "overlapping_with",
        &linear_spheres_type::overlapping_with,
        arg( "object" )
        )
      .def(
        "prefiltered_overlapping_with",
        &linear_spheres_type::template prefiltered_overlapping_with<
          asa::IndexFilter< Vector >
          >,
        ( arg( "object" ), arg( "prefilter" ) )
        )
      .def( "__len__", &linear_spheres_type::size )
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

template < typename Vector >
struct containment_wrappers
{
  static void wrap()
  {
    using namespace boost::python;
    object containment_module(
      handle<>( borrowed( PyImport_AddModule( "asa.containment" ) ) )
      );
    scope().attr( "containment" ) = containment_module;
    scope containment_scope = containment_module;

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
      typename indexing::OverlapFilter< Sphere, overlap::BetweenSpheres >
      >
      filter_helper;

    typedef indexing::PrefilterHelper<
      typename std::vector< Sphere >::const_iterator,
      asa::IndexFilter< vector_type >,
      typename indexing::OverlapFilter< Sphere, overlap::BetweenSpheres >
      >
      prefilter_helper;

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
        "add_from_range",
        add_neighbours_from_range<
          checker_type,
          typename prefilter_helper::filter_range_type >,
        arg( "neighbours" )
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
  containment_wrappers< scitbx::vec3< double > >::wrap();
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

