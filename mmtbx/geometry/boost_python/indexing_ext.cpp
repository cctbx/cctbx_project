#include <boost/python/module.hpp>
#include <boost/python/extract.hpp>
#include <boost/python/class.hpp>

#include <boost/mpl/pair.hpp>
#include <boost/mpl/string.hpp>
#include <boost/mpl/vector.hpp>

#include <boost_adaptbx/exporting.hpp>

#include <scitbx/vec3.h>

#include <mmtbx/geometry/boost_python/indexing.hpp>

namespace mmtbx
{
namespace geometry
{
namespace indexing
{
namespace
{

struct python_code_predicate
{
  python_code_predicate(boost::python::object callable) : m_callable( callable )
  {}

  ~python_code_predicate()
  {}

  bool operator ()(boost::python::object const& object)
  {
    return boost::python::extract< bool >( m_callable( object ) );
  }

private:
  boost::python::object m_callable;
};

template< typename Vector, typename Discrete >
struct python_exports
{
  #pragma clang diagnostic push
  #pragma clang diagnostic ignored "-Wmultichar"
  typedef boost::mpl::vector<
    boost::mpl::pair<
      indexing::Linear< boost::python::object, Vector >,
      boost::mpl::string< 'line', 'ar' >
      >,
    boost::mpl::pair<
      indexing::Hash< boost::python::object, Vector, Discrete >,
      boost::mpl::string< 'hash' >
      >
    > indexers;
  #pragma clang diagnostic pop

  static void wrap()
  {
    boost_adaptbx::exporting::class_list< indexers >::process(
      indexing::python::indexer_exports()
      );
    boost_adaptbx::exporting::class_list< indexers >::process(
      indexing::python::filter_and_range_export< python_code_predicate >(
        "predicate_filtered_"
        )
      );

    using namespace boost::python;
    class_< python_code_predicate >( "predicate", no_init )
      .def( init< boost::python::object >( arg( "callable" ) ) )
      .def( "__call__", &python_code_predicate::operator (), arg( "object" ) )
      ;
  }
};

} // namespace <anonymous>
} // namespace indexing
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_indexing_ext)
{
  mmtbx::geometry::indexing::python_exports< scitbx::vec3< double >, int>::wrap();
}
