/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/eltbx/wavelengths.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace eltbx { namespace wavelengths {
namespace boost_python {

namespace {

  struct characteristic_wrappers
  {
    typedef characteristic w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("characteristic", no_init)
        .def(init<std::string const&>())
        .def("label", &w_t::label)
        .def("as_angstrom", &w_t::as_angstrom)
        .def("as_kev", &w_t::as_kev)
        .def("as_ev", &w_t::as_ev)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    characteristic_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      characteristic,
      characteristic_iterator>::wrap(
        "characteristic_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::wavelengths::boost_python

BOOST_PYTHON_MODULE(wavelengths_ext)
{
  cctbx::eltbx::wavelengths::boost_python::init_module();
}
