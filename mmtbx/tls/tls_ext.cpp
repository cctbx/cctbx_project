#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/tls/tls.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace mmtbx { namespace tls {
namespace {
  boost::python::tuple
  getinitargs(params<> const& self)
  {
    return boost::python::make_tuple(self.t, self.l, self.s);
  }

  void init_module()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_;
   // def("uaniso_from_tls",uaniso_from_tls)
   //;
    class_<uaniso_from_tls>("uaniso_from_tls",
                             init<sym_mat3<double> const&,
                                  sym_mat3<double> const&,
                                  mat3<double> const&,
                                  vec3<double> const&,
                                  vec3<double> const&>())
      .def("u", &uaniso_from_tls::u)
    ;
    class_<tls_from_uaniso_target_and_grads>("tls_from_uaniso_target_and_grads",
                             init<sym_mat3<double> const&,
                                  sym_mat3<double> const&,
                                  mat3<double> const&,
                                  vec3<double> const&,
                                  af::shared<vec3<double> > const&,
                                  af::shared<sym_mat3<double> > const&>())
      .def("target", &tls_from_uaniso_target_and_grads::target)
      .def("grad_T", &tls_from_uaniso_target_and_grads::grad_T)
      .def("grad_L", &tls_from_uaniso_target_and_grads::grad_L)
      .def("grad_S", &tls_from_uaniso_target_and_grads::grad_S)
    ;
    class_<tls_xray_target_grads>("tls_xray_target_grads",
                             init<af::const_ref<cctbx::miller::index<> > const&,
                                  af::const_ref<double> const&,
                                  af::const_ref< std::complex<double> > const&,
                                  cctbx::uctbx::unit_cell const&,
                                  af::const_ref<vec3<double> > const&,
                                  af::shared<sym_mat3<double> > const&,
                                  cctbx::xray::scattering_type_registry const&,
                                  af::const_ref<cctbx::xray::scatterer<> > const&,
                                  af::const_ref<int> const& >())
      .def("target",&tls_xray_target_grads::target)
      .def("gradTLS", &tls_xray_target_grads::gradTLS)
    ;
    typedef return_value_policy<return_by_value> rbv;
    class_<params<> >("params")
      .def(init<scitbx::sym_mat3<double> const&,
                scitbx::sym_mat3<double> const&,
                scitbx::mat3<double> const& >())
      .add_property("t", make_getter(&params<>::t, rbv()))
      .add_property("l", make_getter(&params<>::l, rbv()))
      .add_property("s", make_getter(&params<>::s, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)
    ;
    typedef return_value_policy<return_by_value> rbv;
    class_<tlso<> >("tlso")
      .def(init<scitbx::sym_mat3<double> const&,
                scitbx::sym_mat3<double> const&,
                scitbx::mat3<double> const&,
                scitbx::vec3<double> const& >((arg_("t"),arg_("l"),arg_("s"),
                                               arg_("origin"))))
      .add_property("t",      make_getter(&tlso<>::t,      rbv()))
      .add_property("l",      make_getter(&tlso<>::l,      rbv()))
      .add_property("s",      make_getter(&tlso<>::s,      rbv()))
      .add_property("origin", make_getter(&tlso<>::origin, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)
    ;
    scitbx::af::boost_python::shared_wrapper<params<> >::wrap("shared_params");
    class_<gT_gL_gS>("gT_gL_gS",
                     init<vec3<double> const&,
                          cctbx::uctbx::unit_cell const& ,
                          vec3<double> const& ,
                          cctbx::miller::index<> const&>())
      .def("grad_T", &gT_gL_gS::grad_T)
      .def("grad_L", &gT_gL_gS::grad_L)
      .def("grad_S", &gT_gL_gS::grad_S)
    ;
    def("uaniso_from_tls_one_group",
         (af::shared<sym_mat3<double> >(*)
               (tlso<double>,
                af::shared<vec3<double> > const&)) uaniso_from_tls_one_group,
                                                          (arg_("tlso"),
                                                           arg_("sites_cart")))
   ;
  }

} // namespace <anonymous>
}} // namespace mmtbx::tls

BOOST_PYTHON_MODULE(mmtbx_tls_ext)
{
  mmtbx::tls::init_module();
}
