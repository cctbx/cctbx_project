#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <iotbx/mtz/column.h>

namespace iotbx { namespace mtz {
namespace {

  template <typename DataType>
  struct data_group_wrappers
  {
    typedef data_group<DataType> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
      ;
    }
  };

  struct observations_group_wrappers
  {
    typedef observations_group w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<data_group<double> > >("observations_group", no_init)
        .add_property("sigmas", make_getter(&w_t::sigmas, rbv()))
      ;
    }
  };

  struct complex_group_wrappers
  {
    typedef complex_group w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("complex_group", no_init)
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
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
        .def("extract_integers", &w_t::extract_integers, (
          (arg_("column_label"))))
        .def("extract_reals", &w_t::extract_reals, (
          (arg_("column_label"))))
        .def("extract_reals_anomalous", &w_t::extract_reals_anomalous, (
          arg_("column_label_plus"), arg_("column_label_minus")))
        .def("extract_hls", &w_t::extract_hls, (
          arg_("column_label_a"),
          arg_("column_label_b"),
          arg_("column_label_c"),
          arg_("column_label_d")))
        .def("extract_observations", &w_t::extract_observations, (
          arg_("column_label_data"),
          arg_("column_label_sigmas")))
        .def("extract_observations_anomalous",
          &w_t::extract_observations_anomalous, (
            arg_("column_label_data_plus"),
            arg_("column_label_sigmas_plus"),
            arg_("column_label_data_minus"),
            arg_("column_label_sigmas_minus")))
        .def("extract_delta_anomalous", &w_t::extract_delta_anomalous, (
          arg_("column_label_f_data"),
          arg_("column_label_f_sigmas"),
          arg_("column_label_d_data"),
          arg_("column_label_d_sigmas")))
        .def("extract_complex", &w_t::extract_complex, (
          arg_("column_label_ampl"),
          arg_("column_label_phi")))
        .def("extract_complex_anomalous", &w_t::extract_complex_anomalous, (
          arg_("column_label_ampl_plus"),
          arg_("column_label_phi_plus"),
          arg_("column_label_ampl_minus"),
          arg_("column_label_phi_minus")))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    def("ccp4_liberr_verbosity", ccp4_liberr_verbosity, (arg_("level")));
    data_group_wrappers<int>::wrap("integer_group");
    data_group_wrappers<double>::wrap("real_group");
    data_group_wrappers<cctbx::hendrickson_lattman<> >::wrap("hl_group");
    observations_group_wrappers::wrap();
    complex_group_wrappers::wrap();
    object_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_object() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
