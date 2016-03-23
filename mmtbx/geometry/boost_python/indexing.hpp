#ifndef MMTBX_GEOMETRY_INDEXING_PYTHON_H
#define MMTBX_GEOMETRY_INDEXING_PYTHON_H

#include <string>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>

#include <boost/range/adaptor/filtered.hpp>
#include <boost_adaptbx/boost_range_python.hpp>

#include <mmtbx/geometry/indexing.hpp>

namespace mmtbx
{

namespace geometry
{

namespace indexing
{

namespace python
{

template< typename Indexer >
struct indexer_specific_exports;

template< typename Object, typename Vector >
struct indexer_specific_exports< Linear< Object, Vector > >
{
  typedef Linear< Object, Vector > indexer_type;
  typedef boost::python::class_< indexer_type > python_class_type;

  static void process(python_class_type& myclass)
  {
    using namespace boost::python;
    myclass.def( init<>() )
      ;
  }
};

template< typename Object, typename Vector, typename Discrete >
struct indexer_specific_exports< Hash< Object, Vector, Discrete > >
{
  typedef Hash< Object, Vector, Discrete > indexer_type;
  typedef boost::python::class_< indexer_type > python_class_type;

  static void process(python_class_type& myclass)
  {
    using namespace boost::python;
    typedef typename indexer_type::voxelizer_type voxelizer_type;
    typedef typename indexer_type::discrete_type discrete_type;
    myclass.def(
      init< const voxelizer_type&, const discrete_type& >(
        ( arg( "voxelizer" ), arg( "margin" ) )
        )
      )
      .def( "cubes", &indexer_type::cubes )
      ;
  }
};

struct indexer_exports
{
  template< typename Export >
  void operator()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first indexer_type;
    typedef typename Export::second name_type;

    using namespace boost::python;
    std::string prefix = std::string( boost::mpl::c_str< name_type >::value );

    boost_adaptbx::python::generic_range_wrapper<
      typename indexer_type::range_type
      >::wrap( ( prefix + "_close_objects_range" ).c_str() );

    boost::python::class_< indexer_type > myindexer( prefix.c_str(), no_init );
    myindexer.def( "add", &indexer_type::add, ( arg( "object" ), arg( "position" ) ) )
      .def(
        "close_to",
        &indexer_type::close_to,
        with_custodian_and_ward_postcall< 0, 1 >(),
        arg( "centre" )
        )
      .def( "__len__", &indexer_type::size )
      ;
    indexer_specific_exports< indexer_type >::process( myindexer );
  }
};

template< typename Predicate >
class filter_and_range_export
{
public:
  typedef Predicate predicate_type;

private:
  std::string predicate_id_;

public:
  filter_and_range_export(const std::string& predicate_id)
    : predicate_id_( predicate_id )
  {};

  ~filter_and_range_export() {};

  template< typename Export >
  void operator ()(boost::mpl::identity< Export > myexport) const
  {
    typedef typename Export::first indexer_type;
    typedef typename indexer_type::range_type range_type;
    typedef boost::filtered_range< predicate_type, range_type > filtered_range_type;
    typedef typename Export::second name_type;

    std::string indexer_id = std::string( boost::mpl::c_str< name_type >::value );

    boost_adaptbx::python::generic_range_wrapper< filtered_range_type >
      ::wrap( ( predicate_id_ + indexer_id + "_close_objects_range" ).c_str() );

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

} // namespace python
} // namespace indexing
} // namespace geometry
} // namespace mmtbx

#endif // MMTBX_GEOMETRY_INDEXING_PYTHON_H
