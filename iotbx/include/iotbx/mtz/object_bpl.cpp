#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <iotbx/mtz/column.h>

namespace iotbx { namespace mtz {
namespace {

  struct observation_arrays_wrappers
  {
    typedef observation_arrays w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("observation_arrays", no_init)
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("sigmas", make_getter(&w_t::sigmas, rbv()))
      ;
    }
  };

  struct object_wrappers
  {
    typedef object w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("object", no_init)
        .def(init<>())
        .def(init<af::const_ref<int> const&>(
          (arg_("n_datasets_for_each_crystal"))))
        .def(init<const char*>((arg_("file_name"))))
        .def("title", &w_t::title)
        .def("history", &w_t::history)
        .def("space_group_name", &w_t::space_group_name)
        .def("point_group_name", &w_t::point_group_name)
        .def("space_group", &w_t::space_group)
        .def("n_batches", &w_t::n_batches)
        .def("n_reflections", &w_t::n_reflections)
        .def("space_group_number", &w_t::space_group_number)
        .def("max_min_resolution", &w_t::max_min_resolution)
        .def("n_crystals", &w_t::n_crystals)
        .def("n_active_crystals", &w_t::n_active_crystals)
        .def("crystals", &w_t::crystals)
        .def("lookup_column", &w_t::lookup_column, (arg_("label")))
        .def("valid_indices", &w_t::valid_indices, (arg_("column_label")))
        .def("valid_values", &w_t::valid_values, (arg_("column_label")))
        .def("valid_integers", &w_t::valid_integers, (arg_("column_label")))
        .def("valid_indices_anomalous", &w_t::valid_indices_anomalous, (
          arg_("column_label_plus"), arg_("column_label_minus")))
        .def("valid_values_anomalous", &w_t::valid_values_anomalous, (
          arg_("column_label_plus"), arg_("column_label_minus")))
        .def("valid_complex", &w_t::valid_complex, (
          arg_("column_label_ampl"), arg_("column_label_phi")))
        .def("valid_complex_anomalous", &w_t::valid_complex_anomalous, (
          arg_("column_label_ampl_plus"),
          arg_("column_label_phi_plus"),
          arg_("column_label_ampl_minus"),
          arg_("column_label_phi_minus")))
        .def("valid_delta_anomalous", &w_t::valid_delta_anomalous, (
          arg_("column_label_f_data"),
          arg_("column_label_f_sigmas"),
          arg_("column_label_d_data"),
          arg_("column_label_d_sigmas")))
        .def("valid_hl", &w_t::valid_hl, (
          arg_("column_label_a"),
          arg_("column_label_b"),
          arg_("column_label_c"),
          arg_("column_label_d")))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    def("ccp4_liberr_verbosity", ccp4_liberr_verbosity, (arg_("level")));
    observation_arrays_wrappers::wrap();
    object_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_object() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
