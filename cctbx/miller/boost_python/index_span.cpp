/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/index_span.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct index_span_wrappers
  {
    typedef index_span w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("index_span", no_init)
        .def(init<af::const_ref<index<> > const&>())
        .def("min", &w_t::min)
        .def("max", &w_t::max)
        .def("abs_range", &w_t::abs_range)
        .def("map_grid", &w_t::map_grid)
        .def("is_in_domain", &w_t::is_in_domain)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_index_span()
  {
    index_span_wrappers::wrap();
  }

}}} // namespace cctbx::miller::boost_python
