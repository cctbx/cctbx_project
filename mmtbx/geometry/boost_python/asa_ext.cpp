#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/import.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/object.hpp>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/string.hpp>

#include <string>

#include <scitbx/vec3.h>
#include <boost_adaptbx/boost_range_python.hpp>
#include <boost_adaptbx/exporting.hpp>

#include <mmtbx/geometry/asa.hpp>
#include <mmtbx/geometry/indexing.hpp>
#include <mmtbx/geometry/overlap.hpp>
#include <mmtbx/geometry/containment.hpp>
#include <mmtbx/geometry/sphere_surface_sampling.hpp>
#include <mmtbx/geometry/calculator.hpp>

#include <mmtbx/geometry/boost_python/indexing.hpp>
#include <mmtbx/geometry/boost_python/containment.hpp>

namespace mmtbx
{
namespace geometry
{
namespace asa
{

namespace python
{

template< typename Value >
struct transformed_value_extract
{
public:
  typedef Value value_type;

private:
  boost::python::object transformation_;

public:
  transformed_value_extract(boost::python::object transformation)
  : transformation_( transformation )
  {};

  value_type operator () (boost::python::object obj) const
  {
    return boost::python::extract< value_type >( transformation_( obj ) );
  }
};

template< typename Value >
struct array_adaptor : private transformed_value_extract< Value >
{
public:
  typedef Value value_type;
  typedef std::size_t size_type;

  typedef transformed_value_extract< value_type > transformer_type;

private:
  boost::python::object array_;

public:
  array_adaptor(boost::python::object array, boost::python::object transformation)
  : transformer_type( transformation ), array_( array )
  {};

  value_type operator [](size_type const& index) const
  {
    return transformer_type::operator ()( array_[ index ] );
  }

  size_type size() const
  {
    return boost::python::len( array_ );
  }
};

} // namespace python

namespace
{

template< typename Predicate >
struct filtered_range_type
{
    template< typename Export >
    struct apply
    {
      typedef typename Export::first indexer_type;
      typedef typename indexer_type::range_type range_type;
      typedef boost::filtered_range< Predicate, range_type > type;
    };
};

template< typename Traits >
struct python_exports
{
  static void wrap()
  {
    using namespace boost::python;
    typedef typename Traits::sphere_type sphere_type;
    typedef typename Traits::sphere_bases_type sphere_bases_type;
    typedef typename Traits::vector_type vector_type;
    typedef typename Traits::value_type value_type;
    typedef typename Traits::discrete_type discrete_type;

    // base module
    class_< sphere_type, sphere_bases_type >( "sphere", no_init )
      .def(
        init< const vector_type&, const value_type&, const size_t& >(
          ( arg( "centre" ), arg( "radius" ), arg( "index" ) )
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
      ;

    // indexing module
    object indexing_module(
      handle<>( borrowed( PyImport_AddModule( "asa.indexing" ) ) )
      );
    scope().attr( "indexing" ) = indexing_module;


    { // enter indexing namespace
      scope indexing_scope = indexing_module;

      boost_adaptbx::exporting::class_list< typename Traits::indexers >::process(
        indexing::python::indexer_exports()
        );
    } // exit indexing namespace

    // accessibility module
    object accessibility_module(
      handle<>( borrowed( PyImport_AddModule( "asa.accessibility" ) ) )
      );
    scope().attr( "accessibility" ) = accessibility_module;

    { // enter accessibility namespace
      scope accessibility_scope = accessibility_module;

      typedef typename sphere_surface_sampling::GoldenSpiral< vector_type >
        ::storage_type points_range;
      typedef asa::Transform< vector_type > transformation_type;
      typedef boost::transformed_range< transformation_type, points_range >
        transformed_points_range;

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

      typedef asa::OverlapEqualityFilter< sphere_type, overlap::BetweenSpheres >
          predicate_type;

      class_< predicate_type >( "overlap_equality_predicate", no_init )
        .def( init< const sphere_type& >( arg( "object" ) ) )
        .def( "__call__", &predicate_type::operator (), arg( "other" ) )
        ;

      boost_adaptbx::exporting::class_list< typename Traits::indexers >::process(
        indexing::python::filter_and_range_export< predicate_type >(
          "overlap_equality_filtered_"
          )
        );

      typedef typename boost::mpl::transform<
        typename Traits::indexers,
        filtered_range_type< predicate_type >
        >::type filtered_ranges_type;
      boost_adaptbx::exporting::class_list< typename Traits::checkers >::process(
        containment::python::checker_export<
          filtered_ranges_type,
          transformed_points_range
          >()
        );
    } // exit accessibility namespace
  }
};

template< typename Vector, typename Discrete = int >
struct python_export_traits
{
  typedef Sphere< Vector > sphere_type;
  typedef boost::python::bases< primitive::Sphere< Vector > > sphere_bases_type;

  typedef typename sphere_type::vector_type vector_type;
  typedef typename vector_type::value_type value_type;
  typedef Discrete discrete_type;

  // Exported indexers - indexing module
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmultichar"
  typedef boost::mpl::vector<
    boost::mpl::pair<
      indexing::Linear< sphere_type, vector_type >,
      boost::mpl::string< 'line', 'ar', '_sph', 'eres' >
      >,
    boost::mpl::pair<
      indexing::Hash< sphere_type, vector_type, discrete_type >,
      boost::mpl::string< 'hash', '_sph', 'eres' >
      >
    > indexers;

  // Exported checkers - accessibility module
  typedef boost::mpl::vector<
    boost::mpl::pair<
      containment::Checker< sphere_type, containment::PurePythagorean< false > >,
      boost::mpl::string< 'pyth', 'agor', 'ean' >
      >
    > checkers;
  #pragma clang diagnostic pop
};

struct python_calculator_export
{
  static void wrap()
  {
    using namespace boost::python;
    typedef python::array_adaptor< scitbx::vec3< double > > coordinate_adaptor_type;
    typedef python::array_adaptor< double > radius_adaptor_type;

    class_< coordinate_adaptor_type >( "_coordinate_adaptor", no_init )
      .def( init< object, object >( ( arg( "array" ), arg( "transformation" ) ) ) )
      .def( "__getitem__", &coordinate_adaptor_type::operator [] )
      .def( "__len__", &coordinate_adaptor_type::size )
      ;
    class_< radius_adaptor_type >( "_radius_adaptor", no_init )
      .def( init< object, object >( ( arg( "array" ), arg( "transformation" ) ) ) )
      .def( "__getitem__", &radius_adaptor_type::operator [] )
      .def( "__len__", &radius_adaptor_type::size )
      ;

    typedef calculator::SimpleCalculator<
      coordinate_adaptor_type,
      radius_adaptor_type
      > calculator_type;

    class_< calculator_type >( "calculator", no_init )
      .def(
        init< coordinate_adaptor_type, radius_adaptor_type, double, std::size_t, double, int >(
          (
            arg( "coordinate_adaptor" ),
            arg( "radius_adaptor" ),
            arg( "probe" ) = 1.4,
            arg( "sampling_point_count" ) = 960,
            arg( "cubesize" ) = 7.0,
            arg( "margin" ) = 1
            )
          )
        )
      .def( "accessible_surface_points", &calculator_type::accessible_points, arg( "index" ) )
      .def( "accessible_surface_area", &calculator_type::accessible_surface_area, arg( "index" ) )
      .def(
        "is_overlapping_sphere",
        &calculator_type::is_overlapping_sphere,
        ( arg( "centre" ), arg( "radius" ) )
        )
      ;

    typedef calculator::ConstRadiusCalculator< coordinate_adaptor_type >
      const_radius_calculator_type;

    class_< const_radius_calculator_type >( "const_radius_calculator", no_init )
      .def(
        init< coordinate_adaptor_type, double, double, std::size_t, double, int >(
          (
            arg( "coordinate_adaptor" ),
            arg( "radius" ),
            arg( "probe" ) = 1.4,
            arg( "sampling_point_count" ) = 960,
            arg( "cubesize" ) = 7.0,
            arg( "margin" ) = 1
            )
          )
        )
      .def(
        "accessible_surface_points",
        &const_radius_calculator_type::accessible_points,
        arg( "index" )
        )
      .def(
        "accessible_surface_area",
        &const_radius_calculator_type::accessible_surface_area,
        arg( "index" )
        )
      .def(
        "is_overlapping_sphere",
        &const_radius_calculator_type::is_overlapping_sphere,
        ( arg( "centre" ), arg( "radius" ) )
        )
      ;
  }

};

} // namespace <anonymous>
} // namespace asa
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_asa_ext)
{
  // Dependency
  boost::python::object primitives = boost::python::import( "mmtbx_geometry_primitive_ext" );

  mmtbx::geometry::asa::python_exports<
    mmtbx::geometry::asa::python_export_traits< scitbx::vec3< double > >
    >::wrap();
  mmtbx::geometry::asa::python_calculator_export::wrap();
}

