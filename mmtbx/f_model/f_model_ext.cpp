#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/f_model/f_model.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace mmtbx { namespace f_model {
namespace {

  boost::python::tuple
  getinitargs(core<> const& self)
  {
    return boost::python::make_tuple(self.f_calc,
                                     self.f_mask,
                                     self.b_cart,
                                     self.k_sol,
                                     self.b_sol,
                                     self.overall_scale,
                                     self.hkl,
                                     self.uc,
                                     self.f_model,
                                     self.f_bulk,
                                     self.fb_cart,
                                     self.ss);
  }

  void init_module()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_;

    typedef return_value_policy<return_by_value> rbv;
    class_<core<> >("core")
      .def(init<
           af::shared<std::complex<double> >      const&,
           af::shared<std::complex<double> >      const&,
           scitbx::sym_mat3<double>               const&,
           double                                 const&,
           double                                 const&,
           double                                 const&,
           af::const_ref<cctbx::miller::index<> > const&,
           cctbx::uctbx::unit_cell                const&,
           af::shared<double>                     const& >(
                                                        (arg_("f_calc"),
                                                         arg_("f_mask"),
                                                         arg_("b_cart"),
                                                         arg_("k_sol"),
                                                         arg_("b_sol"),
                                                         arg_("overall_scale"),
                                                         arg_("hkl"),
                                                         arg_("uc"),
                                                         arg_("ss"))))
      .add_property("f_calc",        make_getter(&core<>::f_calc,       rbv()))
      .add_property("f_mask",        make_getter(&core<>::f_mask,       rbv()))
      .add_property("b_cart",        make_getter(&core<>::b_cart,       rbv()))
      .add_property("k_sol",         make_getter(&core<>::k_sol,        rbv()))
      .add_property("b_sol",         make_getter(&core<>::b_sol,        rbv()))
      .add_property("overall_scale", make_getter(&core<>::overall_scale,rbv()))
      .add_property("hkl",           make_getter(&core<>::hkl,          rbv()))
      .add_property("uc",            make_getter(&core<>::uc,           rbv()))
      .add_property("f_model",       make_getter(&core<>::f_model,      rbv()))
      .add_property("f_bulk",        make_getter(&core<>::f_bulk,       rbv()))
      .add_property("fb_cart",       make_getter(&core<>::fb_cart,      rbv()))
      .add_property("ss",            make_getter(&core<>::ss,           rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)
    ;
  }

} // namespace <anonymous>
}} // namespace mmtbx::f_model

BOOST_PYTHON_MODULE(mmtbx_f_model_ext)
{
  mmtbx::f_model::init_module();
}
