#include <cctbx/boost_python/flex_fwd.h>

#include <scitbx/boost_python/utils.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/maptbx/bcr/bcr.h>

#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>

namespace cctbx { namespace maptbx { namespace boost_python {

  boost::python::tuple
  getinitargs_bcr_model(bcr_model<> const& self)
  {
    return boost::python::make_tuple(
      self.scatterer,
      self.B,
      self.C,
      self.R);
  }

boost::python::tuple
  getinitargs_bcr_scatterer(bcr_scatterer<> const& self)
  {
    return boost::python::make_tuple(
      self.site_cart,
      self.u_iso,
      self.occ,
      self.radius,
      self.resolution,
      self.mu,
      self.kappa,
      self.nu,
      self.musq,
      self.kappi);
  }

namespace {

  void init_module()
  {
    using namespace boost::python;

   {
      typedef return_value_policy<return_by_value> rbv;
      class_<bcr_scatterer<> >("bcr_scatterer")
        .def(init<cctbx::cartesian<> const&,
                  double const&,
                  double const&,
                  double const&,
                  double const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  af::shared<double> const& >((arg("site_cart"),
                                               arg("u_iso"),
                                               arg("occ"),
                                               arg("radius"),
                                               arg("resolution"),
                                               arg("mu"),
                                               arg("kappa"),
                                               arg("nu"),
                                               arg("musq"),
                                               arg("kappi"))))
        .add_property("site_cart",  make_getter(&bcr_scatterer<>::site_cart,  rbv()))
        .add_property("u_iso",      make_getter(&bcr_scatterer<>::u_iso,      rbv()))
        .add_property("occ",        make_getter(&bcr_scatterer<>::occ,        rbv()))
        .add_property("radius",     make_getter(&bcr_scatterer<>::radius,     rbv()))
        .add_property("resolution", make_getter(&bcr_scatterer<>::resolution, rbv()))
        .add_property("mu",         make_getter(&bcr_scatterer<>::mu,         rbv()))
        .add_property("kappa",      make_getter(&bcr_scatterer<>::kappa,      rbv()))
        .add_property("nu",         make_getter(&bcr_scatterer<>::nu,         rbv()))
        .add_property("musq",       make_getter(&bcr_scatterer<>::musq,       rbv()))
        .add_property("kappi",      make_getter(&bcr_scatterer<>::kappi,      rbv()))
        .enable_pickling()
        .def("__getinitargs__", getinitargs_bcr_scatterer)
      ;
    }



    {
      typedef return_value_policy<return_by_value> rbv;
      class_<bcr_model<> >("bcr_model")
        .def(init<cctbx::xray::scatterer<> const&,
                  af::shared<double> const&,
                  af::shared<double> const&,
                  af::shared<double> const& >((arg("scatterer"),arg("B"),arg("C"),
                                                 arg("R"))))
        .add_property("scatterer", make_getter(&bcr_model<>::scatterer, rbv()))
        .add_property("B",         make_getter(&bcr_model<>::B,         rbv()))
        .add_property("C",         make_getter(&bcr_model<>::C,         rbv()))
        .add_property("R",         make_getter(&bcr_model<>::R,         rbv()))
        .enable_pickling()
        .def("__getinitargs__", getinitargs_bcr_model)
        //.def("atom_radius",   &bcr_model<>::atom_radius)
        //.def("rho",           &bcr_model<>::rho)
      ;
    }

    {
      class_<image<> >("BCRimage", no_init) // why no_init ?
        .def(init<
          cctbx::uctbx::unit_cell const&,
          boost::python::list const&,
          int const& >(
                    (arg("unit_cell"),
                     arg("bcr_models"),
                     arg("step"))))
        .def("fsc",   &image<>::cc)
        .def("d",     &image<>::d)
        .def("d_inv", &image<>::d_inv)
      ;
    }

    {
      class_<calculator<> >("calculator")
        .def(init<bcr_model<double> const& >((arg("bcr_model") )))
        .def("atom_radius",   &calculator<>::atom_radius)
        .def("rho",           &calculator<>::rho)
      ;
    }

  }

} // namespace <anonymous>

}}} // namespace cctbx::maptbx::boost_python

BOOST_PYTHON_MODULE(cctbx_maptbx_bcr_bcr_ext)
{
  cctbx::maptbx::boost_python::init_module();
}
