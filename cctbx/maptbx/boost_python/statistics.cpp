/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2002 Oct: Created (rwgk)
 */

#include <cctbx/maptbx/statistics.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct statistics_wrappers
  {
    typedef statistics<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("statistics", no_init)
        .def(init<af::const_ref<float, af::flex_grid<> > const&>())
        .def(init<af::const_ref<double, af::flex_grid<> > const&>())
        .def("min", &w_t::min)
        .def("max", &w_t::max)
        .def("mean", &w_t::mean)
        .def("mean_sq", &w_t::mean_sq)
        .def("sigma", &w_t::sigma)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_statistics()
  {
    statistics_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
