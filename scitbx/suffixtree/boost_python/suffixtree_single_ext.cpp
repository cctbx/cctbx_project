#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/with_custodian_and_ward.hpp>
#include <boost/python/stl_iterator.hpp>
#include <boost/python/to_python_converter.hpp>

#include <boost/unordered_map.hpp>

#include <boost_adaptbx/boost_range_python.hpp>
#include <scitbx/boost_python/std_pair.h>

#include <scitbx/suffixtree/word.hpp>
#include <scitbx/suffixtree/boost_python/edge_exports.hpp>
#include <scitbx/suffixtree/boost_python/tree_exports.hpp>
#include <scitbx/suffixtree/boost_python/object_extensions.hpp>
#include <scitbx/suffixtree/matching_statistics.hpp>
#include <scitbx/suffixtree/boost_python/iterator_wrapper.hpp>

namespace scitbx
{
namespace suffixtree
{
namespace
{

template< typename Key, typename Value >
class BoostHashMapAdapter
{
public:
  typedef boost::unordered_map< Key, Value > type;
};

template< typename Glyph >
struct python_exports
{
  typedef Glyph glyph_type;
  typedef word::Single< glyph_type > word_type;
  typedef typename word_type::substring_type substring_type;
  typedef typename word_type::index_type index_type;

  typedef std::size_t suffix_label_type;
  typedef Tree< word_type, suffix_label_type, BoostHashMapAdapter > tree_type;
  typedef MSI< tree_type, boost::python::stl_input_iterator< glyph_type > > msi_type;
  typedef scitbx::suffixtree::python::python_iterator< msi_type > msi_python_iterator;

  static msi_python_iterator get_matching_statistics_range(
    tree_type const& tree,
    boost::python::object const& iterable
    )
  {
    typedef boost::python::stl_input_iterator< glyph_type > iter_type;

    iter_type begin( iterable );
    iter_type end;

    return msi_python_iterator(
      msi_type( tree, begin, end),
      msi_type( tree, end, end )
      );
  }

#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"

  static glyph_type getitem(word_type const& word, std::size_t index)
  {
    if ( index < word.size() )
    {
      return word[ index ];
    }
    else
    {
      PyErr_SetString(PyExc_IndexError, "word index out of range");
      boost::python::throw_error_already_set();
    }
  }

#pragma clang diagnostic pop

  static void wrap()
  {
    using namespace boost::python;

    boost_adaptbx::python::generic_range_wrapper< substring_type >
      ::wrap( "substring" );

    const glyph_type& (word_type::*pindexing)(const index_type&) const =
      &word_type::operator [];

    class_< word_type >( "word", no_init )
      .def( init<>() )
      .def( "append", &word_type::push_back, arg( "glyph" ) )
      .def( "length_descriptor", &word_type::length_ptr )
      .def(
        "substring",
        &word_type::substring,
        with_custodian_and_ward_postcall< 0, 1 >(),
        ( arg( "begin" ), arg( "end" ) )
        )
      .def( "__getitem__", getitem, arg( "index" ) )
      .def( "__len__", &word_type::size )
      ;
    python::edge_exports<
      glyph_type,
      index_type,
      typename word_type::const_length_ptr_type,
      suffix_label_type,
      BoostHashMapAdapter
      >::wrap();
    python::tree_exports< word_type, suffix_label_type, BoostHashMapAdapter >::wrap();
    python::ukkonen_builder_exports< tree_type >::wrap();

/*  class_< msi_range >(
      "matching_statistics_iterator",
      no_init
      )
      .def( "__iter__", boost::python::iterator< msi_range >() )
      ; */

    msi_python_iterator::wrap( "matching_statistics_iterator" );

    to_python_converter<
      typename msi_type::position_type,
      scitbx::boost_python::PairToTupleConverter< typename msi_type::position_type > >();

    to_python_converter<
      typename msi_type::value_type,
      scitbx::boost_python::PairToTupleConverter< typename msi_type::value_type > >();

    def(
      "matching_statistics",
      get_matching_statistics_range,
      ( arg( "tree" ), arg( "iterable" ) )
      );
  }
};

} // namespace <anonymous>
} // namespace suffixtree
} // namespace scitbx


BOOST_PYTHON_MODULE(scitbx_suffixtree_single_ext)
{
  // single module
  scitbx::suffixtree::python_exports< boost::python::object >::wrap();
}
