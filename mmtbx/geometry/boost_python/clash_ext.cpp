#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/import.hpp>

#include <boost/range/adaptor/transformed.hpp>
#include <boost/range/adaptor/filtered.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/pair.hpp>
#include <boost/mpl/string.hpp>

#include <cctbx/sgtbx/rt_mx.h>
#include <scitbx/vec3.h>
#include <iotbx/pdb/small_str.h>

#include <mmtbx/geometry/clash.hpp>
#include <mmtbx/geometry/indexing.hpp>
#include <mmtbx/geometry/overlap.hpp>

#include <mmtbx/geometry/boost_python/exporting.hpp>
#include <mmtbx/geometry/boost_python/indexing.hpp>


namespace mmtbx
{
namespace geometry
{
namespace clash
{
namespace
{

template< typename Altloc >
struct altloc_wrapper
{
  typedef boost::shared_ptr< AltlocStrategy< Altloc > > altloc_strategy_type;

  static bool is_interacting_with(
    const altloc_strategy_type& left,
    const altloc_strategy_type& right
    )
  {
    return left->is_interacting_with( *right );
  }

  static bool is_interacting_with_alternate(
    const altloc_strategy_type& left,
    const char* identifier
    )
  {
    return left->is_interacting_with_alternate( Altloc( identifier, true ) );
  }

  static altloc_strategy_type create_regular()
  {
    return altloc_strategy_type( new RegularAltlocStrategy< Altloc >() );
  }

  static altloc_strategy_type create_alternate(const char* identifier)
  {
    return altloc_strategy_type(
      new AlternateAltlocStrategy< Altloc >( Altloc( identifier, true ) )
      );
  }

  static void wrap()
  {
    using namespace boost::python;

    class_< altloc_strategy_type >( "altloc_strategy", no_init )
      .def(
        "is_interacting_with",
        is_interacting_with,
        arg( "other" )
        )
      .def(
        "is_interacting_with_alternate",
         is_interacting_with_alternate,
         arg( "identifier" )
         )
      .def( "regular", create_regular )
      .staticmethod( "regular" )
      .def( "alternate", create_alternate, arg( "identifier" ) )
      .staticmethod( "alternate" )
      ;
  }
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
    typedef typename Traits::identifier_type identifier_type;
    typedef typename Traits::altloc_strategy_type altloc_strategy_type;
    typedef typename Traits::symop_type symop_type;

    altloc_wrapper< typename Traits::altloc_type >::wrap();

    class_< sphere_type, sphere_bases_type >( "sphere", no_init )
      .def(
        init<
          const vector_type&,
          const value_type&,
          const identifier_type&,
          const identifier_type&,
          const altloc_strategy_type&,
          const symop_type&
          >(
            (
              arg( "centre" ),
              arg( "radius" ),
              arg( "molecule" ),
              arg( "atom" ),
              arg( "altloc" ),
              arg( "symop" )
              )
          )
        )
      .add_property(
        "molecule",
        make_function(
          &sphere_type::molecule,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "atom",
        make_function(
          &sphere_type::atom,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "altloc",
        make_function(
          &sphere_type::altloc_strategy,
          return_value_policy< copy_const_reference >()
          )
        )
      .add_property(
        "symop",
        make_function(
          &sphere_type::symop,
          return_value_policy< copy_const_reference >()
          )
        )
      ;

    // indexing module
    object indexing_module(
      handle<>( borrowed( PyImport_AddModule( "asa.indexing" ) ) )
      );
    scope().attr( "indexing" ) = indexing_module;


    { // enter indexing namespace
      scope indexing_scope = indexing_module;

      exporting::class_list<
        typename Traits::indexers,
        indexing::python::indexer_exports
        >::process();
    } // exit indexing namespace

    typedef OverlapInteractionFilter<
      sphere_type,
      overlap::BetweenSpheresTolerance
      >
      predicate_type;

    class_< predicate_type >( "overlap_interaction_predicate", no_init )
      .def( init< const sphere_type&, const value_type& >(
        ( arg( "object" ), arg( "tolerance" ) ) )
        )
      .def( "__call__", &predicate_type::operator (), arg( "other" ) )
      ;

    exporting::class_list<
      typename Traits::indexers,
      indexing::python::filter_and_range_export< predicate_type >
      >::process(
          indexing::python::filter_and_range_export< predicate_type >(
            "overlap_interaction_filtered_"
            )
        );
  }
};

template< typename Vector, typename Identifier, typename Altloc, typename SymOp >
struct python_export_traits
{
  typedef Sphere< Vector, Identifier, Altloc, SymOp > sphere_type;
  typedef boost::python::bases< primitive::Sphere< Vector > > sphere_bases_type;

  typedef Altloc altloc_type;
  typedef int discrete_type;

  typedef typename sphere_type::vector_type vector_type;
  typedef typename sphere_type::value_type value_type;
  typedef typename sphere_type::identifier_type identifier_type;
  typedef typename sphere_type::altloc_strategy_type altloc_strategy_type;
  typedef typename sphere_type::symop_type symop_type;

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
  #pragma clang diagnostic pop
};

} // namespace <anonymous>
} // namespace clash
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_clash_ext)
{
  // Dependency
  boost::python::object primitives = boost::python::import( "mmtbx_geometry_primitive_ext" );

  mmtbx::geometry::clash::python_exports<
    mmtbx::geometry::clash::python_export_traits<
      scitbx::vec3< double >,
      unsigned long,
      iotbx::pdb::str1,
      cctbx::sgtbx::rt_mx
      >
    >::wrap();
}
