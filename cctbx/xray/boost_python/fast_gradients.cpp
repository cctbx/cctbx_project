#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/xray/fast_gradients.h>
#include <boost/python/class.hpp>

namespace cctbx { namespace xray { namespace boost_python {

namespace {

  struct fast_gradients_wrappers
  {
    typedef fast_gradients<> w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<w_t::base_t> >("fast_gradients", no_init)
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scatterer<> > const&,
                  af::const_ref<double,
                                maptbx::c_grid_padded_p1<3> > const&,
                  af::const_ref<std::complex<double>,
                                maptbx::c_grid_padded_p1<3> > const&,
                  gradient_flags const&,
                  optional<double const&,
                           double const&,
                           double const&,
                           bool> >())
        .def("d_target_d_site_cart", &w_t::d_target_d_site_cart)
        .def("d_target_d_u_iso", &w_t::d_target_d_u_iso)
        .def("d_target_d_u_cart", &w_t::d_target_d_u_cart)
        .def("d_target_d_occupancy", &w_t::d_target_d_occupancy)
        .def("d_target_d_fp", &w_t::d_target_d_fp)
        .def("d_target_d_fdp", &w_t::d_target_d_fdp)
      ;
    }
  };

} // namespace <anoymous>

  void wrap_fast_gradients()
  {
    fast_gradients_wrappers::wrap();
  }

}}} // namespace cctbx::xray::boost_python
