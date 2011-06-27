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
                                     self.f_mask(),
                                     self.k_sols(),
                                     self.b_sol,
                                     self.f_part1,
                                     self.f_part2,
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

    boost::python::to_python_converter<
      af::small< af::shared<std::complex<double> >, max_n_shells> ,
      scitbx::boost_python::container_conversions::to_tuple<
        af::small< af::shared<std::complex<double> >, max_n_shells>   > >();

    scitbx::boost_python::container_conversions::from_python_sequence<
      af::small< af::shared<std::complex<double> >, max_n_shells>,
      scitbx::boost_python::container_conversions::fixed_capacity_policy>();

    typedef return_value_policy<return_by_value> rbv;
    class_<core<> >("core")
      .def(init<
           af::shared<std::complex<double> >      const&,
           af::small< af::shared<std::complex<double> >, max_n_shells > const&,
           af::small< double, max_n_shells>       const&,
           double                                 const&,
           af::shared<std::complex<double> >      const&,
           af::shared<std::complex<double> >      const&,
           scitbx::sym_mat3<double>               const&,
           af::shared<cctbx::miller::index<> >    const&,
           cctbx::uctbx::unit_cell                const&,
           af::shared<double>                     const& >(
                                                        (arg("f_calc"),
                                                         arg("shell_f_masks"),
                                                         arg("k_sols"),
                                                         arg("b_sol"),
                                                         arg("f_part1"),
                                                         arg("f_part2"),
                                                         arg("u_star"),
                                                         arg("hkl"),
                                                         arg("uc"),
                                                         arg("ss"))))
      .add_property("f_calc",        make_getter(&core<>::f_calc,       rbv()))
      .add_property("b_sol",         make_getter(&core<>::b_sol,        rbv()))
      .add_property("f_part1",       make_getter(&core<>::f_part1,      rbv()))
      .add_property("f_part2",       make_getter(&core<>::f_part2,      rbv()))
      .add_property("u_star",        make_getter(&core<>::u_star,       rbv()))
      .add_property("hkl",           make_getter(&core<>::hkl,          rbv()))
      .add_property("uc",            make_getter(&core<>::uc,           rbv()))
      .add_property("f_model",       make_getter(&core<>::f_model,      rbv()))
      .add_property("f_bulk",        make_getter(&core<>::f_bulk,       rbv()))
      .add_property("f_aniso",       make_getter(&core<>::f_aniso,      rbv()))
      .add_property("f_b_sol",       make_getter(&core<>::f_b_sol,      rbv()))
      .add_property("ss",            make_getter(&core<>::ss,           rbv()))
      .def("n_shells", &core<>::n_shells)
      .def("k_sol", &core<>::k_sol)
      .def("k_sols", &core<>::k_sols)
      .def("f_mask", &core<>::f_mask)
      .def("shell_f_mask", &core<>::shell_f_mask)   // XXX re-name to f_mask
      .def("shell_f_masks", &core<>::shell_f_masks) // XXX re-name to f_masks
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
