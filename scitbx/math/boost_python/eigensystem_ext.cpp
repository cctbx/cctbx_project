/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Mar: Created (R.W. Grosse-Kunstleve)
 */

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/boost_python/utils.h>
#include <scitbx/math/eigensystem.h>
#include <boost/python/module.hpp>
#include <boost/python/scope.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>

namespace scitbx { namespace math { namespace eigensystem {
namespace boost_python { namespace {

  struct real_symmetric_wrappers
  {
    typedef real_symmetric<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("real_symmetric", no_init)
        .def(init<af::const_ref<double, af::c_grid<2> > const&,
                  optional<double> >())
        .def(init<scitbx::sym_mat3<double> const&,
                  optional<double> >())
        .def("vectors", &w_t::vectors)
        .def("values", &w_t::values)
      ;
    }
  };

  void init_module()
  {
    using namespace boost::python;

    scope().attr("__version__") = scitbx::boost_python::cvs_revision(
      "$Revision$");

    real_symmetric_wrappers::wrap();
  }

}}}}} // namespace scitbx::math::eigensystem::boost_python::<anonymous>

BOOST_PYTHON_MODULE(eigensystem_ext)
{
  scitbx::math::eigensystem::boost_python::init_module();
}
