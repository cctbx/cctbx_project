#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/cross_correl_components.h>

#include <boost/python/class.hpp>
#include <boost/python/object.hpp>

namespace cctbx { namespace miller { namespace boost_python {

namespace {

  struct cross_correl_components_wrappers
  {
    typedef cross_correlation_components w_t;
    static void wrap()
    {
      using namespace boost::python;
      class_<w_t>("cross_correl_components", no_init)
        .def(init<af::shared<index<> > const&,
                  af::shared<index<> > const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  xfel::hkl_resolution_bins_cpp::hkl_bin_map_type const&,
                  int>())
        .def("N", &w_t::get_count)
        .def("sum_xx", &w_t::get_sum_xx)
        .def("sum_xy", &w_t::get_sum_xy)
        .def("sum_yy", &w_t::get_sum_yy)
        .def("sum_x", &w_t::get_sum_x)
        .def("sum_y", &w_t::get_sum_y)
        ;


    }
  };

} //namespace <anonymous>

  void wrap_cross_correl_components()
  {
    cross_correl_components_wrappers::wrap();
  }


}}} // namespace cctbx::miller::boost_python
