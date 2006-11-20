#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/asu.h>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct index_table_layout_adaptor_wrappers
  {
    typedef index_table_layout_adaptor w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t,bases<sym_equiv_index> >("index_table_layout_adaptor",no_init)
        .def("h", &w_t::h)
        .def("i_column", &w_t::i_column)
      ;
    }
  };

  struct asym_index_wrappers
  {
    typedef asym_index w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<sym_equiv_index> >("asym_index", no_init)
        .def(init<sgtbx::space_group const&,
                  sgtbx::reciprocal_space::asu const&,
                  index<> const&>())
        .def(init<sgtbx::space_group const&,
                  index<> const&>())
        .def(init<sym_equiv_indices const&>())
        .def("one_column", &w_t::one_column)
        .def("two_column", &w_t::two_column)
      ;
    }
  };

  template <typename ValueType>
  struct map_to_asu_no_bool
  {
    static
    void
    wrapper(
      sgtbx::space_group_type const& sg_type,
      bool anomalous_flag,
      af::ref<index<> > const& miller_indices,
      af::ref<ValueType> const& data)
    {
      map_to_asu(sg_type, anomalous_flag, miller_indices, data);
    }
  };

  template <typename ValueType>
  struct map_to_asu_with_bool
  {
    static
    void
    wrapper(
      sgtbx::space_group_type const& sg_type,
      bool anomalous_flag,
      af::ref<index<> > const& miller_indices,
      af::ref<ValueType> const& data,
      bool deg)
    {
      map_to_asu(sg_type, anomalous_flag, miller_indices, data, deg);
    }
  };

} // namespace <anoymous>

  void wrap_asu()
  {
    using namespace boost::python;
    index_table_layout_adaptor_wrappers::wrap();
    asym_index_wrappers::wrap();
    def("map_to_asu",
      (void(*)(sgtbx::space_group_type const&,
               bool,
               af::ref<index<> > const&))
      map_to_asu);
    def("map_to_asu", map_to_asu_no_bool<double>::wrapper);
    def("map_to_asu", map_to_asu_with_bool<double>::wrapper);
    def("map_to_asu", map_to_asu_no_bool<std::complex<double> >::wrapper);
    def("map_to_asu", map_to_asu_no_bool<hendrickson_lattman<> >::wrapper);
    def("is_unique_set_under_symmetry", is_unique_set_under_symmetry, (
      arg_("space_group_type"),
      arg_("anomalous_flag"),
      arg_("miller_indices")));
    def("unique_under_symmetry_selection", unique_under_symmetry_selection, (
      arg_("space_group_type"),
      arg_("anomalous_flag"),
      arg_("miller_indices")));
  }

}}} // namespace cctbx::miller::boost_python
