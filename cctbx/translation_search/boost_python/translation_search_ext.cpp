/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>

namespace cctbx { namespace translation_search { namespace boost_python {

  void wrap_fast_nv1995();
  void wrap_symmetry_flags();

namespace {

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    wrap_fast_nv1995();
    wrap_symmetry_flags();
  }

} // namespace <anonymous>
}}} // namespace cctbx::translation_search::boost_python

BOOST_PYTHON_MODULE(translation_search_ext)
{
  cctbx::translation_search::boost_python::init_module();
}
