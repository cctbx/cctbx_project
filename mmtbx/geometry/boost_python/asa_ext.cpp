#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/operators.hpp>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/string.hpp>

#include <scitbx/vec3.h>
#include <boost_adaptbx/boost_range_python.hpp>

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

template< typename Vector >
Sphere< Vector >
copy_sphere(const Sphere< Vector >& sphere)
{
  return Sphere< Vector >( sphere );
}

template < typename Vector >
struct asa_wrappers
{
  typedef Sphere< Vector > sphere_type;
  typedef primitive::Sphere< Vector > base_sphere_type;

  static void wrap()
  {
    using namespace boost::python;

    typedef typename primitive::Traits< Vector >::vector_type vector_type;
    typedef typename primitive::Traits< Vector >::value_type value_type;

    class_< sphere_type, bases< base_sphere_type > >(
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
          &sphere_type::index,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property( "low", make_function( &sphere_type::low ) )
      .add_property( "high", make_function( &sphere_type::high ) )
      .def( "__hash__", hash_value< sphere_type > )
      .def( self == self )
      .def( self != self )
      .def( "copy", copy_sphere< Vector >, arg( "sample" ) )
      .staticmethod( "copy" )
      ;
  }
};

} // namespace <anonymous>
} // namespace asa

namespace indexing
{

namespace
{

template< typename Sphere, typename Discrete >
struct indexing_wrappers
{
public:
  typedef Sphere sphere_type;
  typedef typename Sphere::vector_type vector_type;
  typedef Linear< Sphere > linear_spheres_type;
  typedef Hash< Sphere, Discrete> hash_spheres_type;
  typedef typename hash_spheres_type::voxelizer_type voxelizer_type;

  typedef boost::mpl::string< '_', 's', 'p', 'h', 'e', 'r', 'e', 's' > suffix_type;
  typedef boost::mpl::string< 'l', 'i', 'n', 'e', 'a', 'r' > linear_prefix_type;
  typedef boost::mpl::string< 'h', 'a', 's', 'h' > hash_prefix_type;

  typedef boost::mpl::vector<
    boost::mpl::pair<
      linear_spheres_type,
      typename boost::mpl::insert_range<
        linear_prefix_type,
        typename boost::mpl::end< linear_prefix_type >::type,
        suffix_type
        >::type
      >,
    boost::mpl::pair<
      hash_spheres_type,
      typename boost::mpl::insert_range<
        hash_prefix_type,
        typename boost::mpl::end< hash_prefix_type >::type,
        suffix_type
        >::type
      >
    > exports;

  static void wrap()
  {
    using namespace boost::python;
    object indexing_module(
      handle<>( borrowed( PyImport_AddModule( "asa.indexing" ) ) )
      );
    scope().attr( "indexing" ) = indexing_module;
    scope indexing_scope = indexing_module;

    wrap_linear_indexer();
    wrap_hash_indexer();
  }

  template< typename Indexer >
  void
  static wrap_output_range(const std::string& prefix)
  {
    typedef typename Indexer::range_type range_type;
    boost_adaptbx::python::generic_range_wrapper< range_type >::wrap(
      ( prefix + "_close_objects_range" ).c_str()
      );
  }

  static void wrap_linear_indexer()
  {
    wrap_output_range< linear_spheres_type >( "linear_spheres" );

    using namespace boost::python;

    class_< linear_spheres_type >( "linear_spheres", no_init )
      .def( init<>() )
      .def( "add", &linear_spheres_type::add, arg( "object" ) )
      .def(
        "close_to",
        &linear_spheres_type::close_to,
        with_custodian_and_ward_postcall< 0, 1 >(),
        arg( "object" )
        )
      .def( "__len__", &linear_spheres_type::size )
      ;
  }

  static void wrap_hash_indexer()
  {
    wrap_output_range< hash_spheres_type >( "hash_spheres" );

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
        with_custodian_and_ward_postcall< 0, 1 >(),
        arg( "object" )
        )
      .def( "__len__", &hash_spheres_type::size )
      .def( "cubes", &hash_spheres_type::cubes )
      ;
  }
};

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

  template< typename ExportType >
  void operator ()(boost::mpl::identity< ExportType > myexport)
  {
    typedef typename boost::mpl::at_c< ExportType, 2 >::type filtered_range_type;
    myclass_.def(
      "add",
      add_neighbours_from_range< Checker, filtered_range_type >,
      boost::python::arg( "neighbours" )
      );
  }
};

template< typename PredicateType >
struct indexer_filtered_range_export_type
{
  template< typename ExportType >
  struct apply
  {
    typedef typename ExportType::first indexer_type;
    typedef typename indexer_type::range_type range_type;
    typedef boost::filtered_range< PredicateType, range_type > filtered_range_type;
    typedef boost::mpl::vector<
      PredicateType,
      range_type,
      filtered_range_type,
      typename ExportType::second
      > type;
  };
};

struct filtered_range_python_export
{
  template< typename ExportType >
  void operator ()(boost::mpl::identity< ExportType > myrange)
  {
    typedef typename boost::mpl::at_c< ExportType, 0 >::type predicate_type;
    typedef typename boost::mpl::at_c< ExportType, 1 >::type range_type;
    typedef typename boost::mpl::at_c< ExportType, 2 >::type filtered_range_type;
    typedef typename boost::mpl::at_c< ExportType, 3 >::type prefix_type;

    std::string prefix = std::string( boost::mpl::c_str< prefix_type >:: value );

    boost_adaptbx::python::generic_range_wrapper< filtered_range_type >
      ::wrap( ( "filtered_" + prefix + "_close_objects_range" ).c_str() );

    filtered_range_type
      (*filterfunc)( range_type&, predicate_type ) =
        &boost::adaptors::filter< range_type, predicate_type >;

    using namespace boost::python;

    def(
      "filter",
      filterfunc,
      with_custodian_and_ward_postcall< 0, 1 >(),
      ( arg( "range" ), arg( "predicate" ) )
      );
    }
};

template < typename IndexerWrapper >
struct accessibility_wrappers
{
  typedef typename IndexerWrapper::sphere_type sphere_type;
  typedef typename sphere_type::vector_type vector_type;
  typedef typename sphere_surface_sampling::GoldenSpiral< vector_type >
    ::storage_type points_range;
  typedef asa::Transform< vector_type > transformation_type;
  typedef boost::transformed_range< transformation_type, points_range >
    transformed_points_range;
  typedef asa::OverlapEqualityFilter< sphere_type, overlap::BetweenSpheres >
    predicate_type;
  typedef typename boost::mpl::transform<
    typename IndexerWrapper::exports,
    indexer_filtered_range_export_type< predicate_type >
    >::type filtered_range_exports;

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

    class_< predicate_type >( "overlap_equality_predicate", no_init )
      .def( init< const sphere_type& >( arg( "object" ) ) )
      .def( "__call__", &predicate_type::operator (), arg( "other" ) )
      ;

    boost::mpl::for_each<
      filtered_range_exports,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( filtered_range_python_export() );

    typedef Checker< sphere_type, PurePythagorean< false > >
      pythagorean_checker;
    boost_adaptbx::python::generic_range_wrapper< typename pythagorean_checker::storage_type >
      ::wrap( "spheres_range" );

    create_wrappings< pythagorean_checker >( "pythagorean" );
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
      filtered_range_exports,
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

} // namespace <anonymous>
} // namespace containment
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_asa_ext)
{
  typedef mmtbx::geometry::asa::asa_wrappers< scitbx::vec3< double > > asa_wrapper;
  asa_wrapper::wrap();
  typedef mmtbx::geometry::indexing::indexing_wrappers<
    typename asa_wrapper::sphere_type,
    int
    > indexing_wrapper;
  indexing_wrapper::wrap();
  mmtbx::geometry::containment::accessibility_wrappers< indexing_wrapper >::wrap();
}

