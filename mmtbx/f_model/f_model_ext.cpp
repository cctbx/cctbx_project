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
                                     self.k_sol,
                                     self.b_sol,
                                     self.f_part,
                                     self.k_part,
                                     self.b_part,
                                     self.u_star,
                                     self.hkl,
                                     self.uc,
                                     self.f_model,
                                     self.f_bulk,
                                     self.f_aniso,
                                     self.f_b_sol,
                                     self.ss);
  }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;

    typedef return_value_policy<return_by_value> rbv;
    class_<core<> >("core")
      .def(init<
           af::shared<std::complex<double> >      const&,
           af::shared<std::complex<double> >      const&,
           double                                 const&,
           double                                 const&,
           af::shared<std::complex<double> >      const&,
           double                                 const&,
           double                                 const&,
           scitbx::sym_mat3<double>               const&,
           af::shared<cctbx::miller::index<> >    const&,
           cctbx::uctbx::unit_cell                const&,
           af::shared<double>                     const& >(
                                                        (arg("f_calc"),
                                                         arg("f_mask"),
                                                         arg("k_sol"),
                                                         arg("b_sol"),
                                                         arg("f_part_base"),
                                                         arg("k_part"),
                                                         arg("b_part"),
                                                         arg("u_star"),
                                                         arg("hkl"),
                                                         arg("uc"),
                                                         arg("ss"))))
      .add_property("f_calc",        make_getter(&core<>::f_calc,       rbv()))
      .add_property("f_mask",        make_getter(&core<>::f_mask,       rbv()))
      .add_property("k_sol",         make_getter(&core<>::k_sol,        rbv()))
      .add_property("b_sol",         make_getter(&core<>::b_sol,        rbv()))
      .add_property("f_part",        make_getter(&core<>::f_part,       rbv()))
      .add_property("k_part",        make_getter(&core<>::k_part,       rbv()))
      .add_property("b_part",        make_getter(&core<>::b_part,       rbv()))
      .add_property("u_star",        make_getter(&core<>::u_star,       rbv()))
      .add_property("hkl",           make_getter(&core<>::hkl,          rbv()))
      .add_property("uc",            make_getter(&core<>::uc,           rbv()))
      .add_property("f_model",       make_getter(&core<>::f_model,      rbv()))
      .add_property("f_bulk",        make_getter(&core<>::f_bulk,       rbv()))
      .add_property("f_aniso",       make_getter(&core<>::f_aniso,      rbv()))
      .add_property("f_b_sol",       make_getter(&core<>::f_b_sol,      rbv()))
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
