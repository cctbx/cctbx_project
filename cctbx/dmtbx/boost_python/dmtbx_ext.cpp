/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Nov: Created (Ralf W. Grosse-Kunstleve)
 */

#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/class.hpp>

namespace cctbx { namespace dmtbx { namespace boost_python {

  void wrap_triplet();

namespace {

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    scitbx::boost_python::import_module(
      "cctbx_boost.array_family.flex_cctbx_ext");

    wrap_triplet();
  }

} // namespace <anonymous>
}}} // namespace cctbx::dmtbx::boost_python

BOOST_PYTHON_MODULE(dmtbx_ext)
{
  cctbx::dmtbx::boost_python::init_module();
}
