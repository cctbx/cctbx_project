/* Copyright (c) 2001-2002 The Regents of the University of California
   through E.O. Lawrence Berkeley National Laboratory, subject to
   approval by the U.S. Department of Energy.
   See files COPYRIGHT.txt and LICENSE.txt for further details.

   Revision history:
     2003 Jun: Created (rwgk)
 */

#include <scitbx/array_family/boost_python/flex_fwd.h>

#include <scitbx/math/mean_and_variance.h>
#include <boost/python/class.hpp>

namespace scitbx { namespace af { namespace boost_python { namespace {

  struct mean_and_variance_wrappers
  {
    typedef math::mean_and_variance<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("mean_and_variance", no_init)
        .def(init<af::const_ref<double> const&>())
        .def(init<af::const_ref<double> const&,
                  af::const_ref<double> const&>())
        .def("mean", &w_t::mean)
        .def("gsl_variance", &w_t::gsl_variance)
        .def("gsl_standard_deviation", &w_t::gsl_standard_deviation)
        .def("cumulative_variance", &w_t::cumulative_variance)
        .def("cumulative_standard_deviation",
          &w_t::cumulative_standard_deviation)
        .def("conservative_variance", &w_t::conservative_variance)
        .def("conservative_standard_deviation",
          &w_t::conservative_standard_deviation)
        .def("sum_weights", &w_t::sum_weights)
        .def("sum_weights_sq", &w_t::sum_weights_sq)
        .def("sum_weights_values", &w_t::sum_weights_values)
        .def("sum_weights_delta_sq", &w_t::sum_weights_delta_sq)
      ;
    }
  };

} // namespace <anonymous>

  void wrap_flex_mean_and_variance()
  {
    mean_and_variance_wrappers::wrap();
  }

}}} // namespace scitbx::af::boost_python
