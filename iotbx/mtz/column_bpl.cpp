#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_arg.hpp>
#include <iotbx/mtz/column.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>

namespace iotbx { namespace mtz {
namespace {

  struct column_wrappers
  {
    typedef column w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("column", no_init)
        .def(init<dataset const&, int>((
          arg("mtz_dataset"), arg("i_column"))))
        .def("mtz_dataset", &w_t::mtz_dataset)
        .def("i_column", &w_t::i_column)
        .def("mtz_crystal", &w_t::mtz_crystal)
        .def("mtz_object", &w_t::mtz_object)
        .def("label", &w_t::label)
        .def("set_label", &w_t::set_label, (arg("new_label")),
          return_self<>())
        .def("type", &w_t::type)
        .def("set_type", &w_t::set_type, (arg("new_type")), return_self<>())
        .def("is_active", &w_t::is_active)
        .def("source", &w_t::source)
        .def("set_source", &w_t::set_source, (arg("new_source")),
          return_self<>())
        .def("group_name", &w_t::group_name)
        .def("set_group_name", &w_t::set_group_name, (arg("new_group_name")),
          return_self<>())
        .def("group_type", &w_t::group_type)
        .def("set_group_type", &w_t::set_group_type, (arg("new_group_type")),
          return_self<>())
        .def("group_position", &w_t::group_position)
        .def("set_group_position", &w_t::set_group_position, (
          arg("new_group_position")),
            return_self<>())
        .def("array_size", &w_t::array_size)
        .def("array_capacity", &w_t::array_capacity)
        .def("path", &w_t::path)
        .def("get_other", &w_t::get_other, (arg("label")))
        .def("n_valid_values", &w_t::n_valid_values)
        .def("extract_valid_values", &w_t::extract_valid_values)
        .def("selection_valid", &w_t::selection_valid)
        .def("extract_values", &w_t::extract_values, (
          arg("not_a_number_substitute")=0))
        .def("set_values", (void(w_t::*)(
          af::const_ref<float> const&,
          af::const_ref<bool> const&) const)&w_t::set_values, (
            arg("values"), arg("selection_valid")))
        .def("set_values", (void(w_t::*)(
          af::const_ref<float> const&) const)&w_t::set_values, (
            arg("values")))
        .def("set_reals",
          (af::shared<int>(w_t::*)(
            af::const_ref<cctbx::miller::index<> > const&,
            af::const_ref<double> const&))
              &w_t::set_reals, (
          arg("miller_indices"), arg("data")))
        .def("set_reals",
          (void(w_t::*)(
            af::const_ref<int> const&,
            af::const_ref<double> const&))
              &w_t::set_reals, (
          arg("mtz_reflection_indices"), arg("data")))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_column");
      }
    }
  };

  void
  wrap_all()
  {
    column_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_column() { wrap_all(); }

}}} // namespace iotbx::mtz::boost_python
