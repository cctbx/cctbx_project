/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Sep: Created (rwgk)
 */

#include <cctbx/sgtbx/brick.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace sgtbx { namespace boost_python {

namespace {

  struct brick_wrappers
  {
    typedef brick w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("brick", no_init)
        .def(init<space_group_type const&>())
        .def("as_string", &w_t::as_string)
        .def("__str__", &w_t::as_string)
        .def("is_inside", (bool(w_t::*)(tr_vec const&) const) &w_t::is_inside)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_brick()
  {
    brick_wrappers::wrap();
  }

}}} // namespace cctbx::sgtbx::boost_python
