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
#include <boost/version.hpp>
#if BOOST_VERSION >= 103000
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#endif

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
#if BOOST_VERSION >= 103000
      typedef return_value_policy<return_by_value> rbv;
#endif
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
#if BOOST_VERSION >= 103000
        .add_property("label", make_getter(&w_t::label, rbv()),
                               make_setter(&w_t::label))
        .add_property("caasf", make_getter(&w_t::caasf, rbv()),
                               make_setter(&w_t::caasf))
        .add_property("fp_fdp", make_getter(&w_t::fp_fdp, rbv()),
                                make_setter(&w_t::fp_fdp))
        .add_property("site", make_getter(&w_t::site, rbv()),
                              make_setter(&w_t::site))
        .add_property("occupancy", make_getter(&w_t::occupancy, rbv()),
                                   make_setter(&w_t::occupancy))
        .add_property("anisotropic_flag",
          make_getter(&w_t::anisotropic_flag, rbv()),
          make_setter(&w_t::anisotropic_flag))
        .add_property("u_iso", make_getter(&w_t::u_iso, rbv()),
                               make_setter(&w_t::u_iso))
        .add_property("u_star", make_getter(&w_t::u_star, rbv()),
                                make_setter(&w_t::u_star))
#else
        .def_readwrite("label", &w_t::label)
        .def_readwrite("caasf", &w_t::caasf)
        .def_readwrite("fp_fdp", &w_t::fp_fdp)
        .def_readwrite("site", &w_t::site)
        .def_readwrite("occupancy", &w_t::occupancy)
        .def_readwrite("anisotropic_flag", &w_t::anisotropic_flag)
        .def_readwrite("u_iso", &w_t::u_iso)
        .def_readwrite("u_star", &w_t::u_star)
#endif
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
