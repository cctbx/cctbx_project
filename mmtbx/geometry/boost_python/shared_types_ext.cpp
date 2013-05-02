#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/to_python_converter.hpp>

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/identity.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/placeholders.hpp>

#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/algorithm/iteration/for_each.hpp>

#include <boost/bind.hpp>

#include <scitbx/vec3.h>
#include <mmtbx/geometry/indexing.hpp>

namespace mmtbx
{

namespace geometry
{

namespace shared
{

// Adapted from http://mail.python.org/pipermail/cplusplus-sig/attachments/20090227/0dd51fec/attachment.hpp
struct ListBuilder
{
  typedef void result_type;

  template< typename T >
  void operator ()(boost::python::list mylist, const T& value) const
  {
    mylist.append( value );
  }
};

template< typename FusionSequence >
struct FusionSequenceConverter
{
  static PyObject* convert(const FusionSequence& seq)
  {
    boost::python::list mylist;
    boost::fusion::for_each(
      seq,
      boost::bind( ListBuilder(), mylist, _1 )
      );

    return boost::python::incref( boost::python::tuple( mylist ).ptr() );
  }
};

namespace
{

struct to_python_fusion_sequence_export
{
  template< typename ExportType >
  void operator ()(boost::mpl::identity< ExportType > myrange)
  {
    using namespace boost::python;
    to_python_converter<
      ExportType,
      FusionSequenceConverter< ExportType >
      >();
  }
};

struct type_wrappers
{
public:
  typedef boost::mpl::vector<
    boost::fusion::vector3< int, int, int >
    >
    to_python_export_types;

  static void wrap()
  {
    typedef scitbx::vec3< double > vector_type;
    typedef mmtbx::geometry::indexing::Voxelizer<
      vector_type,
      boost::fusion::vector3< int, int, int >,
      int
      > voxelizer_type;

    using namespace boost::python;

    class_< voxelizer_type >( "voxelizer", no_init )
      .def(
        init< const vector_type&, const vector_type& >(
          ( arg( "base" ), arg( "step" ) )
          )
        )
      .def( "__call__", &voxelizer_type::operator (), arg( "vector" ) )
      ;

    boost::mpl::for_each<
      to_python_export_types,
      boost::mpl::make_identity< boost::mpl::placeholders::_ >
      >( to_python_fusion_sequence_export() );
  }
};

} // namespace <anonymous>
} // namespace shared
} // namespace geometry
} // namespace mmtbx

BOOST_PYTHON_MODULE(mmtbx_geometry_shared_types_ext)
{
  mmtbx::geometry::shared::type_wrappers::wrap();
}
