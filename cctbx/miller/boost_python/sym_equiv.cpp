/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/miller/sym_equiv.h>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct sym_equiv_index_wrappers
  {
    typedef sym_equiv_index w_t;
    typedef hendrickson_lattman<> h_l;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      ht_angle_overloads, ht_angle, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      mate_overloads, mate, 0, 1)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      phase_eq_overloads, phase_eq, 1, 2)

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      phase_in_overloads, phase_in, 1, 2)

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
        .def("ht_angle", &w_t::ht_angle, ht_angle_overloads())
        .def("friedel_flag", &w_t::friedel_flag)
        .def("mate", &w_t::mate, mate_overloads())
        .def("phase_eq",
          (double(w_t::*)(double const&, bool) const) 0,
          phase_eq_overloads())
        .def("phase_in",
          (double(w_t::*)(double, bool) const) 0,
          phase_in_overloads())
        .def("complex_eq",
          (std::complex<double>(w_t::*)(std::complex<double> const&) const)
          &w_t::complex_eq)
        .def("complex_in",
          (std::complex<double>(w_t::*)(std::complex<double>) const)
          &w_t::complex_in)
        .def("hendrickson_lattman_eq",
          (h_l(w_t::*) (h_l const&) const) &w_t::hendrickson_lattman_eq)
        .def("hendrickson_lattman_in",
          (h_l(w_t::*) (h_l) const) &w_t::hendrickson_lattman_in)
      ;
    }
  };

  struct sym_equiv_indices_wrappers
  {
    typedef sym_equiv_indices w_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      is_valid_phase_overloads, is_valid_phase, 1, 3)

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("sym_equiv_indices", no_init)
        .def(init<sgtbx::space_group const&, index<> const&>())
        .def("phase_restriction", &w_t::phase_restriction)
        .def("is_centric", &w_t::is_centric)
        .def("indices", &w_t::indices, ccr())
        .def("multiplicity", &w_t::multiplicity)
        .def("f_mates", &w_t::f_mates)
        .def("epsilon", &w_t::epsilon)
        .def("__call__",
          (sym_equiv_index(w_t::*)(std::size_t) const) &w_t::operator())
        .def("is_valid_phase",
          &w_t::is_valid_phase, is_valid_phase_overloads())
        .def("p1_listing", &w_t::p1_listing)
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
