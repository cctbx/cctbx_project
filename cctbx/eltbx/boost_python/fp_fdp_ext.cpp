/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <scitbx/boost_python/utils.h>
#include <cctbx/eltbx/fp_fdp.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace eltbx { namespace boost_python {

namespace {

  struct fp_fdp_wrappers
  {
    typedef fp_fdp w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("fp_fdp", no_init)
        .def(init<optional<float, float> >())
        .def("is_valid_fp", &w_t::is_valid_fp)
        .def("is_valid_fdp", &w_t::is_valid_fdp)
        .def("is_valid", &w_t::is_valid)
        .def("fp", &w_t::fp)
        .def("fdp", &w_t::fdp)
        .def("as_complex", &w_t::as_complex_double)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    fp_fdp_wrappers::wrap();
  }

} // namespace <anonymous>
}}} // namespace cctbx::eltbx::boost_python

BOOST_PYTHON_MODULE(fp_fdp_ext)
{
  cctbx::eltbx::boost_python::init_module();
}
