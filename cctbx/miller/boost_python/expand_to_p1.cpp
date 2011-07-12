#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/expand_to_p1.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct expand_to_p1_iselection_wrappers
  {
    typedef expand_to_p1_iselection w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<index<> > const&,
          bool>((
            arg("space_group"),
            arg("anomalous_flag"),
            arg("indices"),
            arg("build_iselection"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("iselection", make_getter(&w_t::iselection, rbv()))
      ;
    }
  };

  template <typename DataType, typename WrappedType>
  struct expand_to_p1_data_wrappers
  {
    typedef WrappedType w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<index<> > const&,
          af::const_ref<DataType> const&>((
            arg("space_group"),
            arg("anomalous_flag"),
            arg("indices"),
            arg("data"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
      ;
    }
  };

  template <typename FloatType>
  struct expand_to_p1_phases_wrappers
  {
    typedef expand_to_p1_phases<FloatType> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<index<> > const&,
          af::const_ref<FloatType> const&,
          bool>((
            arg("space_group"),
            arg("anomalous_flag"),
            arg("indices"),
            arg("data"),
            arg("deg"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_expand_to_p1()
  {
    using namespace boost::python;
    expand_to_p1_iselection_wrappers::wrap("expand_to_p1_iselection");
    expand_to_p1_data_wrappers<
      std::complex<double>,
      expand_to_p1_complex<double> >::wrap(
        "expand_to_p1_complex");
    expand_to_p1_data_wrappers<
      hendrickson_lattman<double>,
      expand_to_p1_hendrickson_lattman<double> >::wrap(
        "expand_to_p1_hendrickson_lattman");
    expand_to_p1_phases_wrappers<double>::wrap("expand_to_p1_phases");
  }

}}} // namespace cctbx::miller::boost_python
