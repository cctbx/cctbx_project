#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python/return_arg.hpp>
#include <iotbx/mtz/batch.h>
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
        .def_readonly("anomalous_flag", &w_t::anomalous_flag)
        .add_property("mtz_reflection_indices",
          make_getter(&w_t::mtz_reflection_indices, rbv()))
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
        .def_readonly("anomalous_flag", &w_t::anomalous_flag)
        .add_property("mtz_reflection_indices",
          make_getter(&w_t::mtz_reflection_indices, rbv()))
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
        .def(init<const char*>((arg("file_name"))))
        .def("title", &w_t::title)
        .def("set_title", &w_t::set_title, (
          arg("title"), arg("append")=false), return_self<>())
        .def("history", &w_t::history)
        .def("add_history",
          (object&(w_t::*)(af::const_ref<std::string> const&))
            &w_t::add_history, (
          arg("lines")), return_self<>())
        .def("add_history",
          (object&(w_t::*)(const char*))
            &w_t::add_history, (
          arg("line")), return_self<>())
        .def("space_group_name", &w_t::space_group_name)
        .def("set_space_group_name", &w_t::set_space_group_name, (
          arg("name")), return_self<>())
        .def("space_group_number", &w_t::space_group_number)
        .def("set_space_group_number", &w_t::set_space_group_number, (
          arg("number")), return_self<>())
        .def("point_group_name", &w_t::point_group_name)
        .def("set_point_group_name", &w_t::set_point_group_name, (
          arg("name")), return_self<>())
        .def("lattice_centring_type",
          &w_t::lattice_centring_type)
        .def("set_lattice_centring_type",
          &w_t::set_lattice_centring_type, (
            arg("symbol")), return_self<>())
        .def("n_symmetry_matrices", &w_t::n_symmetry_matrices)
        .def("space_group_confidence", &w_t::space_group_confidence)
        .def("space_group", &w_t::space_group)
        .def("set_space_group", &w_t::set_space_group, (
          arg("space_group")), return_self<>())
        .def("reserve", &w_t::reserve)
        .def("n_batches", &w_t::n_batches)
        .def("batches", &w_t::batches)
        .def("add_batch", &w_t::add_batch)
        .def("sort_batches", &w_t::sort_batches)
        .def("n_reflections", &w_t::n_reflections)
        .def("max_min_resolution", &w_t::max_min_resolution)
        .def("n_crystals", &w_t::n_crystals)
        .def("n_active_crystals", &w_t::n_active_crystals)
        .def("crystals", &w_t::crystals)
        .def("add_crystal",
          (crystal(w_t::*)(
            const char*, const char*, af::double6 const&))
              &w_t::add_crystal, (
          arg("name"),
          arg("project_name"),
          arg("unit_cell_parameters")))
        .def("add_crystal",
          (crystal(w_t::*)(
            const char*, const char*, cctbx::uctbx::unit_cell const&))
              &w_t::add_crystal, (
          arg("name"),
          arg("project_name"),
          arg("unit_cell")))
        .def("has_crystal", &w_t::has_crystal, (arg("name")))
        .def("has_column", &w_t::has_column, (arg("label")))
        .def("get_column", &w_t::get_column, (arg("label")))
        .def("extract_miller_indices", &w_t::extract_miller_indices)
        .def("replace_miller_indices", &w_t::replace_miller_indices, (
          (arg("miller_indices"))))
        .def("extract_integers",
          (integer_group(w_t::*)(const char*) const) &w_t::extract_integers, (
          (arg("column_label"))))
        .def("extract_integers",
          (af::shared<int>(w_t::*)
            (af::const_ref<int> const&, const char*) const)
              &w_t::extract_integers, (
          (arg("mtz_reflection_indices"), arg("column_label"))))
        .def("extract_integers_anomalous", &w_t::extract_integers_anomalous, (
          arg("column_label_plus"), arg("column_label_minus")))
        .def("extract_reals",
          (real_group(w_t::*)(const char*) const) &w_t::extract_reals, (
          (arg("column_label"))))
        .def("extract_reals",
          (af::shared<double>(w_t::*)
            (af::const_ref<int> const&, const char*) const)
              &w_t::extract_reals, (
          (arg("mtz_reflection_indices"), arg("column_label"))))
        .def("extract_reals_anomalous", &w_t::extract_reals_anomalous, (
          arg("column_label_plus"), arg("column_label_minus")))
        .def("extract_hendrickson_lattman",
            &w_t::extract_hendrickson_lattman, (
          arg("column_label_a"),
          arg("column_label_b"),
          arg("column_label_c"),
          arg("column_label_d")))
        .def("extract_hendrickson_lattman_ab_only",
            &w_t::extract_hendrickson_lattman_ab_only, (
          arg("column_label_a"),
          arg("column_label_b")))
        .def("extract_hendrickson_lattman_anomalous",
            &w_t::extract_hendrickson_lattman_anomalous, (
          arg("column_label_a_plus"),
          arg("column_label_b_plus"),
          arg("column_label_c_plus"),
          arg("column_label_d_plus"),
          arg("column_label_a_minus"),
          arg("column_label_b_minus"),
          arg("column_label_c_minus"),
          arg("column_label_d_minus")))
        .def("extract_hendrickson_lattman_anomalous_ab_only",
            &w_t::extract_hendrickson_lattman_anomalous_ab_only, (
          arg("column_label_a_plus"),
          arg("column_label_b_plus"),
          arg("column_label_a_minus"),
          arg("column_label_b_minus")))
        .def("extract_observations", &w_t::extract_observations, (
          arg("column_label_data"),
          arg("column_label_sigmas")))
        .def("extract_observations_anomalous",
          &w_t::extract_observations_anomalous, (
            arg("column_label_data_plus"),
            arg("column_label_sigmas_plus"),
            arg("column_label_data_minus"),
            arg("column_label_sigmas_minus")))
        .def("extract_delta_anomalous", &w_t::extract_delta_anomalous, (
          arg("column_label_f_data"),
          arg("column_label_f_sigmas"),
          arg("column_label_d_data"),
          arg("column_label_d_sigmas")))
        .def("extract_complex", &w_t::extract_complex, (
          arg("column_label_ampl"),
          arg("column_label_phi")))
        .def("extract_complex_anomalous", &w_t::extract_complex_anomalous, (
          arg("column_label_ampl_plus"),
          arg("column_label_phi_plus"),
          arg("column_label_ampl_minus"),
          arg("column_label_phi_minus")))
        .def("write", &w_t::write, (arg("file_name")))
        .def("xml", &w_t::xml)
        .def("unknown_headers", &w_t::unknown_headers)
        .def("number_of_unknown_headers", &w_t::number_of_unknown_headers)
        .def("delete_reflection", &w_t::delete_reflection, (arg("iref")))
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    def("cmtz_struct_sizes", cmtz_struct_sizes);
    def("ccp4_liberr_verbosity", ccp4_liberr_verbosity, (arg("level")));
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
