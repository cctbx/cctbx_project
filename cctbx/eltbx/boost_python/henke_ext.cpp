/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/eltbx/henke.h>
#include <scitbx/boost_python/utils.h>
#include <scitbx/boost_python/iterator_wrappers.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace eltbx { namespace henke { namespace boost_python {

namespace {

  struct table_wrappers
  {
    typedef table w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("table", no_init)
        .def(init<std::string const&, optional<bool> >())
        .def("label", &w_t::label)
        .def("atomic_number", &w_t::atomic_number)
        .def("at_ev", &w_t::at_ev)
        .def("at_kev", &w_t::at_kev)
        .def("at_angstrom", &w_t::at_angstrom)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    scitbx::boost_python::import_module(
      "cctbx_boost.eltbx.fp_fdp_ext");

    table_wrappers::wrap();
    scitbx::boost_python::iterator_wrappers<
      table, table_iterator>::wrap("table_iterator");
  }

} // namespace <anonymous>
}}}} // namespace cctbx::eltbx::henke::boost_python

BOOST_PYTHON_MODULE(henke_ext)
{
  cctbx::eltbx::henke::boost_python::init_module();
}
