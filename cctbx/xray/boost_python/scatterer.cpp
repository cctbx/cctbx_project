/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/scatterer.h>
#include <boost/python/class.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct scatterer_wrappers
  {
    typedef scatterer<> w_t;
    typedef w_t::float_type flt_t;

    BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(
      apply_symmetry_overloads, apply_symmetry, 2, 6)

    static w_t
    copy(w_t const& o) { return o; }

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("scatterer", no_init)
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  flt_t const&,
                  flt_t const&,
                  eltbx::caasf::wk1995 const&,
                  std::complex<flt_t> const&>())
        .def(init<std::string const&,
                  fractional<flt_t> const&,
                  scitbx::sym_mat3<flt_t> const&,
                  flt_t const&,
                  eltbx::caasf::wk1995 const&,
                  std::complex<flt_t> const&>())
        .def_readwrite("label", &w_t::label)
        .def_readwrite("caasf", &w_t::caasf)
        .def_readwrite("fp_fdp", &w_t::fp_fdp)
        .def_readwrite("site", &w_t::site)
        .def_readwrite("occupancy", &w_t::occupancy)
        .def_readwrite("anisotropic_flag", &w_t::anisotropic_flag)
        .def_readwrite("u_iso", &w_t::u_iso)
        .def_readwrite("u_star", &w_t::u_star)
        .def("apply_symmetry", &w_t::apply_symmetry,apply_symmetry_overloads())
        .def("update_weight", &w_t::update_weight)
        .def("multiplicity", &w_t::multiplicity)
        .def("weight", &w_t::weight)
        .def("copy", copy)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_scatterer()
  {
    scatterer_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
