#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/sampled_model_density.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct sampled_model_density_wrappers
  {
    typedef sampled_model_density<> w_t;
    typedef w_t::grid_point_type grid_point_type;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<w_t::base_t> >("sampled_model_density", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scatterer<> > const&,
                  grid_point_type const&,
                  grid_point_type const&,
                  optional<double const&,
                           double const&,
                           double const&,
                           bool,
                           bool> >())
        .def("real_map", &w_t::real_map)
        .def("complex_map", &w_t::complex_map)
        .def("eliminate_u_extra_and_normalize",
          &w_t::eliminate_u_extra_and_normalize)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_sampled_model_density()
  {
    sampled_model_density_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
