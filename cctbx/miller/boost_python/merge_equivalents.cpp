/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Apr: Created (R.W. Grosse-Kunstleve)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/merge_equivalents.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct merge_equivalents_wrappers
  {
    typedef merge_equivalents<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("merge_equivalents", no_init)
        .def(init<af::const_ref<index<> > const&,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&>())
        .def("indices", &w_t::indices)
        .def("data", &w_t::data)
        .def("sigmas", &w_t::sigmas)
        .def("redundancies", &w_t::redundancies)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_merge_equivalents()
  {
    merge_equivalents_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
