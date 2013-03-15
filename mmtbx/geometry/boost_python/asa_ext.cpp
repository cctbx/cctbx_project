#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>
#include <boost/python/list.hpp>
#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>

#include <scitbx/vec3.h>
#include <boost_adaptbx/boost_range_python.hpp>

#include <mmtbx/geometry/asa.hpp>
#include <mmtbx/geometry/indexing.hpp>
#include <mmtbx/geometry/overlap.hpp>
#include <mmtbx/geometry/containment.hpp>
#include <mmtbx/geometry/sphere_surface_sampling.hpp>
#include <mmtbx/geometry/range_adaptors.hpp>

#include <string>

namespace scitbx
{

size_t
hash_value(const scitbx::vec3< int >& voxel)
{
  std::size_t seed = 0;
  boost::hash_combine( seed, voxel[0] );
  boost::hash_combine( seed, voxel[1] );
  boost::hash_combine( seed, voxel[2] );

  return seed;
}

} // namespace scitbx

namespace mmtbx
{
namespace geometry
{
namespace asa
{
namespace
{

template< typename Vector >
Sphere< Vector >
copy_sphere(const Sphere< Vector >& sphere)
{
  return Sphere< Vector >( sphere );
}

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
        init< const vector_type&, const value_type& >(
          ( arg( "centre" ), arg( "radius" ) )
          )
        )
      .add_property(
        "index",
        make_function(
          &Sphere< Vector >::index,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property( "low", make_function( &Sphere< Vector >::low ) )
      .add_property( "high", make_function( &Sphere< Vector >::high ) )
      .def( "__hash__", hash_value< Sphere< Vector > > )
      .def( self == self )
      .def( self != self )
      .def( "copy", copy_sphere< Vector >, arg( "sample" ) )
      .staticmethod( "copy" )
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

template < typename Predicate, typename InputRange >
struct filtered_range_type
{
  typedef boost::filtered_range< Predicate, InputRange > type;
};

template < typename Sphere, typename Voxel >
struct indexing_wrappers
{
public:
  typedef Sphere sphere_type;
  typedef typename Sphere::vector_type vector_type;
  typedef Linear< Sphere > linear_spheres_type;
  typedef Voxelizer< vector_type, Voxel > voxelizer_type;
  typedef Hash< Sphere, voxelizer_type > hash_spheres_type;
  typedef adaptors::OverlapEqualityFilter< Sphere, overlap::BetweenSpheres >
    predicate_type;

  typedef boost::mpl::vector< linear_spheres_type, hash_spheres_type > indexers;

  static void wrap()
  {
    using namespace boost::python;
    object indexing_module(
      handle<>( borrowed( PyImport_AddModule( "asa.indexing" ) ) )
      );
    scope().attr( "indexing" ) = indexing_module;
    scope indexing_scope = indexing_module;

    class_< predicate_type >( "overlap_equality_predicate", no_init )
      .def( init< const Sphere& >( arg( "object" ) ) )
      .def( "__call__", &predicate_type::operator (), arg( "other" ) )
      ;

    typedef typename indexing::Linear< sphere_type >::range_type
      linear_close_objects_range;

    wrap_linear_indexer();
    wrap_hash_indexer();
  }

  template< typename Indexer >
  void
  static wrap_ranges_and_filters(const std::string& prefix)
  {
    typedef typename Indexer::range_type range_type;
    boost_adaptbx::python::generic_range_wrapper< range_type >::wrap(
      ( prefix + "_close_objects_range" ).c_str()
      );

    typedef typename filtered_range_type< predicate_type, range_type >::type
      filtered_range;
    boost_adaptbx::python::generic_range_wrapper< filtered_range >
      ::wrap( ( "filtered_" + prefix + "_close_objects_range" ).c_str() );

    using namespace boost::python;

    boost::filtered_range< predicate_type, range_type >
      (*filterfunc)( range_type&, predicate_type ) =
        &boost::adaptors::filter< range_type, predicate_type >;

    def(
      "filter",
      filterfunc,
      with_custodian_and_ward_postcall< 0, 1 >(),
      ( arg( "range" ), arg( "predicate" ) )
      );
  }

  static void wrap_linear_indexer()
  {
    wrap_ranges_and_filters< linear_spheres_type >( "linear_spheres" );

    using namespace boost::python;

    class_< linear_spheres_type >( "linear_spheres", no_init )
      .def( init<>() )
      .def( "add", &linear_spheres_type::add, arg( "object" ) )
      .def(
        "close_to",
        &linear_spheres_type::close_to,
        return_internal_reference<>(),
        arg( "object" )
        )
      .def( "__len__", &linear_spheres_type::size )
      ;
  }

  static void wrap_hash_indexer()
  {
    wrap_ranges_and_filters< hash_spheres_type >( "hash_spheres" );

    using namespace boost::python;
    class_< voxelizer_type >( "voxelizer", no_init )
      .def(
        init< const vector_type&, const vector_type& >(
          ( arg( "base" ), arg( "step" ) )
          )
        )
      .def( "__call__", &voxelizer_type::operator (), arg( "vector" ) )
      ;

    class_< hash_spheres_type >( "hash_spheres", no_init )
      .def( init< const voxelizer_type& >( arg( "voxelizer" ) ) )
      .def( "add", &hash_spheres_type::add, arg( "object" ) )
      .def(
        "close_to",
        &hash_spheres_type::close_to,
        arg( "object" )
        )
      .def( "__len__", &hash_spheres_type::size )
      .def( "cubes", &hash_spheres_type::cubes )
      .def( "count", &hash_spheres_type::count )
      ;
  }
};

void init_module()
{
  typedef asa::Sphere< scitbx::vec3< double > > sphere_type;
  indexing_wrappers< sphere_type, scitbx::vec3< int > >::wrap();
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

template< typename Checker, typename InputRange >
void
add_neighbours_from_range(Checker& checker, const InputRange& neighbours)
{
  checker.add( neighbours.begin(), neighbours.end() );
}

template< typename Checker >
class add_method_definer
{
public:
  typedef boost::python::class_< Checker > python_class;

private:
  python_class myclass_;

public:
  add_method_definer(const python_class& myclass) : myclass_( myclass ) {};
  ~add_method_definer() {};

  template< typename RangeType >
  void operator ()(boost::mpl::identity< RangeType > myrangetype)
  {
    myclass_.def(
      "add",
      add_neighbours_from_range< Checker, RangeType >,
      boost::python::arg( "neighbours" )
      );
  }
};

template< typename Predicate >
struct indexer_filtered_range_type
{
  template< typename Indexer >
  struct apply
  {
    typedef boost::filtered_range< Predicate, typename Indexer::range_type >
      type;
  };
};

template < typename IndexerWrapper >
struct accessibility_wrappers
{
  typedef typename IndexerWrapper::sphere_type sphere_type;
  typedef typename sphere_type::vector_type vector_type;
  typedef typename sphere_surface_sampling::GoldenSpiral< vector_type >
    ::storage_type points_range;
  typedef adaptors::Transform< vector_type > transformation_type;
  typedef boost::transformed_range< transformation_type, points_range >
    transformed_points_range;
  typedef typename boost::mpl::transform<
    typename IndexerWrapper::indexers,
    indexer_filtered_range_type< typename IndexerWrapper::predicate_type >
    >::type filtered_range_types;

  static void wrap()
  {
    using namespace boost::python;
    object accessibility_module(
      handle<>( borrowed( PyImport_AddModule( "asa.accessibility" ) ) )
      );
    scope().attr( "accessibility" ) = accessibility_module;
    scope accessibility_scope = accessibility_module;

    boost_adaptbx::python::generic_range_wrapper< transformed_points_range >
      ::wrap( "transformed_points_range" );

    class_< transformation_type >( "transformation", no_init )
      .def(
        init< const vector_type&, const typename vector_type::value_type& >(
          ( arg( "centre" ), arg( "radius" ) )
          )
        )
      .def( "__call__", &transformation_type::operator (), arg( "point" ) )
      ;

    boost::transformed_range< transformation_type, points_range >
      (*transformfunc)( points_range&, transformation_type ) =
        &boost::adaptors::transform< transformation_type, points_range >;

    def(
      "transform",
      transformfunc,
      with_custodian_and_ward_postcall< 0, 1 >(),
      ( arg( "range" ), arg( "transformation" ) )
      );

    typedef Checker< sphere_type, PurePythagorean< false > >
      pythagorean_checker;
    create_wrappings< pythagorean_checker >( "pythagorean" );

    /*
    create_wrappings< sphere_type, DiamondPrefilter< false > >(
      "diamond_prefilter"
      );
    */
  }

  template< typename Checker >
  static void create_wrappings(const std::string& name)
  {
    using namespace boost::python;

    class_< Checker > myclass( ( name + "_checker" ).c_str(), no_init );
    myclass \
      .def( init<>() )
      .def(
        "add",
        add_neighbours_from_list< Checker >,
        arg( "neighbours" )
        )
      .def(
        "neighbours",
        &Checker::neighbours,
        return_internal_reference<>()
        )
      .def(
        "__call__",
        &Checker::operator (),
        arg( "point" )
        )
      ;

    boost::mpl::for_each<
      filtered_range_types,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( add_method_definer< Checker >( myclass ) );

    typedef boost::filtered_range< Checker, transformed_points_range >
      filtered_transformed_points_range;
    boost_adaptbx::python::generic_range_wrapper<
      filtered_transformed_points_range
      >
      ::wrap( "filtered_transformed_points_range" );

    boost::filtered_range< Checker, transformed_points_range >
      (*filterfunc)( transformed_points_range&, Checker ) =
        &boost::adaptors::filter< transformed_points_range, Checker >;

    def(
      "filter",
      filterfunc,
      with_custodian_and_ward_postcall< 0, 1 >(),
      ( arg( "range" ), arg( "predicate" ) )
      );
  }
};

void init_module()
{
  typedef asa::Sphere< scitbx::vec3< double > > sphere_type;
  typedef scitbx::vec3< int > voxel_type;
  typedef indexing::indexing_wrappers< sphere_type, voxel_type >
    indexer_wrapper_type;
  accessibility_wrappers< indexer_wrapper_type >::wrap();
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

