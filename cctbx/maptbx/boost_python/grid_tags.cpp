#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/grid_tags.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace maptbx {
namespace {

  struct grid_tags_wrappers
  {
    typedef grid_tags<> w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      dependent_correlation_overloads, dependent_correlation, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      verify_overloads, verify, 1, 2)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      class_<w_t>("grid_tags", no_init)
        .def(init<af::int3 const&>((arg_("dim"))))
        .def("is_valid", &w_t::is_valid)
        .def("tag_array", &w_t::tag_array)
        .def("build", &w_t::build,
          (arg_("space_group_type"), arg_("symmetry_flags")))
        .def("space_group_type", &w_t::space_group_type, rir())
        .def("symmetry_flags", &w_t::symmetry_flags, ccr())
        .def("grid_ss", &w_t::grid_ss, ccr())
        .def("n_grid_misses", &w_t::n_grid_misses)
        .def("n_independent", &w_t::n_independent)
        .def("n_dependent", &w_t::n_dependent)
        .def("dependent_correlation",
          (scitbx::math::linear_correlation<>(w_t::*)(
            af::const_ref<float, af::c_grid_padded<3> > const&,
            double) const)
              &w_t::dependent_correlation,
              dependent_correlation_overloads(
                (arg_("data"), arg_("epsilon"))))
        .def("dependent_correlation",
          (scitbx::math::linear_correlation<>(w_t::*)(
            af::const_ref<double, af::c_grid_padded<3> > const&,
            double) const)
              &w_t::dependent_correlation,
              dependent_correlation_overloads(
                (arg_("data"), arg_("epsilon"))))
        .def("verify",
          (bool(w_t::*)(
            af::const_ref<float, af::c_grid_padded<3> > const&,
            double) const)
              &w_t::verify, verify_overloads(
                (arg_("data"), arg_("min_correlation"))))
        .def("verify",
          (bool(w_t::*)(
            af::const_ref<double, af::c_grid_padded<3> > const&,
            double) const)
              &w_t::verify, verify_overloads(
                (arg_("data"), arg_("min_correlation"))))
        .def("sum_sym_equiv_points",
          (void(w_t::*)(af::ref<float, c_grid_padded_p1<3> > const&) const)
            &w_t::sum_sym_equiv_points,
              (arg_("data")))
        .def("sum_sym_equiv_points",
          (void(w_t::*)(af::ref<double, c_grid_padded_p1<3> > const&) const)
            &w_t::sum_sym_equiv_points,
              (arg_("data")))
      ;
    }
  };

} // namespace <anoymous>

namespace boost_python {

  void wrap_grid_tags()
  {
    grid_tags_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
