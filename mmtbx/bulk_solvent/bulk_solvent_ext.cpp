#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/bulk_solvent/bulk_solvent.h>

namespace mmtbx { namespace bulk_solvent {
namespace {

  void init_module()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_;
    class_<target_gradients_aniso>("target_gradients_aniso",
                             init<af::const_ref<double> const&,
                               af::const_ref< std::complex<double> > const&,
                               af::const_ref< std::complex<double> > const&,
                               sym_mat3<double> const&,
                               double const&,
                               double const&,
                               af::const_ref<cctbx::miller::index<> > const&,
                               cctbx::uctbx::unit_cell const&,
                               bool const&,
                               bool const&,
                               bool const&>())
      .def("target", &target_gradients_aniso::target)
      .def("grad_b_cart", &target_gradients_aniso::grad_b_cart)
      .def("scale_target", &target_gradients_aniso::scale_target)
      .def("grad_ksol", &target_gradients_aniso::grad_ksol)
      .def("grad_bsol", &target_gradients_aniso::grad_bsol)
    ;
    class_<target_gradients_aniso_ml>("target_gradients_aniso_ml",
                             init<af::const_ref<double> const&,
                                  af::const_ref< std::complex<double> > const&,
                                  af::const_ref< std::complex<double> > const&,
                                  sym_mat3<double> const&,
                                  double const&,
                                  double const&,
                                  af::const_ref<cctbx::miller::index<> > const&,
                                  cctbx::uctbx::unit_cell const&,
                                  cctbx::sgtbx::space_group const&,
                                  af::const_ref<bool> const&,
                                  af::const_ref<double> const&,
                                  af::const_ref<double> const&,
                                  double const&>())
      .def("target", &target_gradients_aniso_ml::target)
      .def("grad_b_cart", &target_gradients_aniso_ml::grad_b_cart)
      .def("grad_ksol", &target_gradients_aniso_ml::grad_ksol)
      .def("grad_bsol", &target_gradients_aniso_ml::grad_bsol)
      .def("grad_k", &target_gradients_aniso_ml::grad_k)
    ;
    def("r_factor",r_factor)
   ;
    def("scale",scale)
   ;
    def("fb_cart",fb_cart)
   ;
    def("r_factor_aniso_fast",r_factor_aniso_fast)
   ;
    def("ksol_bsol_grid_search",ksol_bsol_grid_search)
   ;
   def("symmetrize_mask",
      (void(*)
        (af::ref<int, af::c_grid<3> > const&,
         af::const_ref<long, af::c_grid<3> > const&)) symmetrize_mask);
  }

} // namespace <anonymous>
}} // namespace mmtbx::bulk_solvent

BOOST_PYTHON_MODULE(mmtbx_bulk_solvent_ext)
{
  mmtbx::bulk_solvent::init_module();
}
