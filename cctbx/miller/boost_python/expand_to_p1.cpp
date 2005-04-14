#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/expand_to_p1.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct expand_to_p1_indices_wrappers
  {
    typedef expand_to_p1_indices w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("expand_to_p1_indices", no_init)
        .def(init<
          sgtbx::space_group const&,
          bool,
          af::const_ref<index<> > const&>((
            arg_("space_group"),
            arg_("anomalous_flag"),
            arg_("indices"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
      ;
    }
  };

  template <typename ComplexType, typename WrappedType>
  struct expand_to_p1_generic_wrappers
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
          af::const_ref<ComplexType> const&>((
            arg_("space_group"),
            arg_("anomalous_flag"),
            arg_("indices"),
            arg_("data"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
      ;
    }
  };

  template <typename FloatType>
  struct expand_to_p1_obs_wrappers
  {
    typedef expand_to_p1_obs<FloatType> w_t;

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
          af::const_ref<FloatType> const&>((
            arg_("space_group"),
            arg_("anomalous_flag"),
            arg_("indices"),
            arg_("data"),
            arg_("sigmas"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("sigmas", make_getter(&w_t::sigmas, rbv()))
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
            arg_("space_group"),
            arg_("anomalous_flag"),
            arg_("indices"),
            arg_("data"),
            arg_("deg"))))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_expand_to_p1()
  {
    expand_to_p1_indices_wrappers::wrap();
    expand_to_p1_generic_wrappers<
      bool,
      expand_to_p1_scalar<bool> >::wrap(
        "expand_to_p1_bool");
    expand_to_p1_generic_wrappers<
      int,
      expand_to_p1_scalar<int> >::wrap(
        "expand_to_p1_int");
    expand_to_p1_generic_wrappers<
      double,
      expand_to_p1_scalar<double> >::wrap(
        "expand_to_p1_double");
    expand_to_p1_generic_wrappers<
      std::complex<double>,
      expand_to_p1_complex<double> >::wrap(
        "expand_to_p1_complex");
    expand_to_p1_generic_wrappers<
      hendrickson_lattman<double>,
      expand_to_p1_hendrickson_lattman<double> >::wrap(
        "expand_to_p1_hendrickson_lattman");
    expand_to_p1_obs_wrappers<double>::wrap("expand_to_p1_obs");
    expand_to_p1_phases_wrappers<double>::wrap("expand_to_p1_phases");
  }

}}} // namespace cctbx::miller::boost_python
