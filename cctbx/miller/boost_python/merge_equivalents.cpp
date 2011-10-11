#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/hendrickson_lattman.h>
#include <cctbx/miller/merge_equivalents.h>
#include <boost/python/class.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct merge_equivalents_real_wrappers
  {
    typedef merge_equivalents_real<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("merge_equivalents_real", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<double> const&>())
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
        .add_property("r_linear", make_getter(&w_t::r_linear, rbv()))
        .add_property("r_square", make_getter(&w_t::r_square, rbv()))
        .add_property("r_int", &w_t::r_int)
      ;
    }
  };

  struct merge_equivalents_complex_wrappers
  {
    typedef merge_equivalents_generic<std::complex<double>, double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("merge_equivalents_complex", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<std::complex<double> > const&>())
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
      ;
    }
  };

  struct merge_equivalents_hl_wrappers
  {
    typedef merge_equivalents_generic<hendrickson_lattman<>, double> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("merge_equivalents_hl", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<hendrickson_lattman<> > const&>())
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
      ;
    }
  };

  template <typename IntegralType>
  struct merge_equivalents_exact_wrappers
  {
    typedef merge_equivalents_exact<IntegralType> w_t;

    static void
    wrap(const char* python_name)
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>(python_name, no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<IntegralType> const&>())
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
      ;
    }
  };

  struct merge_equivalents_obs_wrappers
  {
    typedef merge_equivalents_obs<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("merge_equivalents_obs", no_init)
        .def(init<
          af::const_ref<index<> > const&,
          af::const_ref<double> const&,
          af::const_ref<double> const&,
          double>((
            arg("unmerged_indices"),
            arg("unmerged_data"),
            arg("unmerged_sigmas"),
            arg("sigma_dynamic_range")=1e-6)))
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("sigmas", make_getter(&w_t::sigmas, rbv()))
        .def_readonly("sigma_dynamic_range", &w_t::sigma_dynamic_range)
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
        .add_property("r_linear", make_getter(&w_t::r_linear, rbv()))
        .add_property("r_square", make_getter(&w_t::r_square, rbv()))
        .add_property("r_int", &w_t::r_int)
      ;
    }
  };

  struct merge_equivalents_shelx_wrappers
  {
    typedef merge_equivalents_shelx<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("merge_equivalents_shelx", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&>())
        .add_property("indices", make_getter(&w_t::indices, rbv()))
        .add_property("data", make_getter(&w_t::data, rbv()))
        .add_property("sigmas", make_getter(&w_t::sigmas, rbv()))
        .add_property("redundancies", make_getter(&w_t::redundancies, rbv()))
        .add_property("r_linear", make_getter(&w_t::r_linear, rbv()))
        .add_property("r_square", make_getter(&w_t::r_square, rbv()))
        .add_property("r_int", &w_t::r_int)
        .add_property("inconsistent_equivalents", &w_t::inconsistent_equivalents)
      ;
    }
  };
} // namespace <anoymous>

  void wrap_merge_equivalents()
  {
    merge_equivalents_real_wrappers::wrap();
    merge_equivalents_complex_wrappers::wrap();
    merge_equivalents_hl_wrappers::wrap();
    merge_equivalents_exact_wrappers<bool>::wrap(
      "merge_equivalents_exact_bool");
    merge_equivalents_exact_wrappers<int>::wrap(
      "merge_equivalents_exact_int");
    merge_equivalents_obs_wrappers::wrap();
    merge_equivalents_shelx_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
