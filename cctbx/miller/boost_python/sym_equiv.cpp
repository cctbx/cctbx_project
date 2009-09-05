#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/sym_equiv.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct sym_equiv_index_wrappers
  {
    typedef sym_equiv_index w_t;
    typedef hendrickson_lattman<> h_l;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("sym_equiv_index", no_init)
        .def("h", &w_t::h)
        .def("hr", &w_t::hr, ccr())
        .def("ht", &w_t::ht)
        .def("t_den", &w_t::t_den)
        .def("ht_angle", &w_t::ht_angle, (arg("deg")=false))
        .def("friedel_flag", &w_t::friedel_flag)
        .def("mate", &w_t::mate, (arg("i_mate")=1))
        .def("phase_eq",
          (double(w_t::*)(double const&, bool) const) &w_t::phase_eq, (
            arg("phi_in"), arg("deg")=false))
        .def("phase_in",
          (double(w_t::*)(double, bool) const) &w_t::phase_in, (
            arg("phi_eq"), arg("deg")=false))
        .def("complex_eq",
          (std::complex<double>(w_t::*)(std::complex<double> const&) const)
          &w_t::complex_eq, (arg("f_in")))
        .def("complex_in",
          (std::complex<double>(w_t::*)(std::complex<double> const&) const)
          &w_t::complex_in, (arg("f_eq")))
        .def("hendrickson_lattman_eq",
          (h_l(w_t::*) (h_l const&) const) &w_t::hendrickson_lattman_eq, (
            arg("hl_in")))
        .def("hendrickson_lattman_in",
          (h_l(w_t::*) (h_l) const) &w_t::hendrickson_lattman_in, (
            arg("hl_eq")))
      ;
    }
  };

  struct sym_equiv_indices_wrappers
  {
    typedef sym_equiv_indices w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("sym_equiv_indices", no_init)
        .def(init<sgtbx::space_group const&, index<> const&>((
          arg("space_group"), arg("h_in"))))
        .def("phase_restriction", &w_t::phase_restriction)
        .def("is_centric", &w_t::is_centric)
        .def("indices", &w_t::indices, ccr())
        .def("multiplicity", &w_t::multiplicity, (arg("anomalous_flag")))
        .def("f_mates", &w_t::f_mates, (arg("anomalous_flag")))
        .def("epsilon", &w_t::epsilon)
        .def("__call__",
          (sym_equiv_index(w_t::*)(std::size_t) const) &w_t::operator(), (
            arg("i")))
        .def("is_valid_phase", &w_t::is_valid_phase, (
          arg("phi"), arg("deg")=false, arg("tolerance")=1e-5))
        .def("p1_listing", &w_t::p1_listing, (arg("anomalous_flag")))
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sym_equiv()
  {
    sym_equiv_index_wrappers::wrap();
    sym_equiv_indices_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
