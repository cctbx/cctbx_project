#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_arg.hpp>
#include <boost/python/overloads.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/stl/map_wrapper.h>
#include <scitbx/stl/vector_wrapper.h>
#include <cctbx/crystal/pair_tables.h>

namespace cctbx { namespace crystal {
namespace {

  struct pair_sym_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<pair_sym_dict, rir>::wrap(
        "pair_sym_dict");
      scitbx::af::boost_python::shared_wrapper<pair_sym_dict, rir>::wrap(
        "pair_sym_table");
    }
  };

  struct pair_asu_table_table_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_internal_reference<> rir;
      scitbx::stl::boost_python::map_wrapper<pair_asu_dict, rir>::wrap(
        "pair_asu_dict");
      scitbx::af::boost_python::shared_wrapper<pair_asu_dict, rir>::wrap(
        "pair_asu_table_table");
    }
  };

  struct pair_asu_table_wrappers
  {
    typedef pair_asu_table<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      add_all_pairs_overloads, add_all_pairs, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      extract_pair_sym_table_overloads, extract_pair_sym_table, 0, 1)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("pair_asu_table", no_init)
        .def(init<
          boost::shared_ptr<direct_space_asu::asu_mappings<> >&>(
            (arg_("asu_mappings"))))
        .def("asu_mappings", &w_t::asu_mappings)
        .def("table", &w_t::table, ccr())
        .def("__contains__",
          (bool(w_t::*)(direct_space_asu::asu_mapping_index_pair const&))
            &w_t::contains, (
          arg_("pair")))
        .def("contains",
          (bool(w_t::*)(unsigned, unsigned, unsigned))
            &w_t::contains, (
          arg_("i_seq"), arg_("j_seq"), arg_("j_sym")))
        .def("__eq__", &w_t::operator==)
        .def("__ne__", &w_t::operator!=)
        .def("add_all_pairs", &w_t::add_all_pairs,
          add_all_pairs_overloads((
            arg_("distance_cutoff"), arg_("epsilon")=1.e-6))[return_self<>()])
        .def("add_pair_sym_table", &w_t::add_pair_sym_table, (
          arg_("sym_table")), return_self<>())
        .def("add_pair", (pair_asu_table<>&(w_t::*)(
            unsigned, unsigned, sgtbx::rt_mx const&)) &w_t::add_pair,
          (arg_("i_seq"), arg_("j_seq"), arg_("rt_mx_ji")), return_self<>())
        .def("add_pair", (pair_asu_table<>&(w_t::*)(
            af::tiny<unsigned, 2> const&)) &w_t::add_pair,
          (arg_("i_seqs")), return_self<>())
        .def("extract_pair_sym_table", &w_t::extract_pair_sym_table,
          extract_pair_sym_table_overloads((
            arg_("skip_j_seq_less_than_i_seq")=true)))
      ;
    }
  };

  void
  wrap_all()
  {
    pair_sym_table_wrappers::wrap();
    pair_asu_table_table_wrappers::wrap();
    pair_asu_table_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_pair_tables() { wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
