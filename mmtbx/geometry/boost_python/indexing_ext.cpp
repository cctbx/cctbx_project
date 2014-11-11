#include <boost/python/module.hpp>

#include <boost/mpl/pair.hpp>
#include <boost/mpl/string.hpp>

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

template< typename Vector, typename Discrete >
struct python_exports
{
  static void wrap()
  {
    #pragma clang diagnostic push
    #pragma clang diagnostic ignored "-Wmultichar"
    python::indexer_exports indexer_exporter;
    indexer_exporter(
      boost::mpl::identity<
        boost::mpl::pair<
          Linear< boost::python::object, Vector >,
          boost::mpl::string< 'line', 'ar' >
          >
        >()
      );
    indexer_exporter(
      boost::mpl::identity<
        boost::mpl::pair<
          Hash< boost::python::object, Vector, Discrete >,
          boost::mpl::string< 'hash' >
          >
        >()
      );
    #pragma clang diagnostic pop
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
