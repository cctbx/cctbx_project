#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/rotamer/fit.h>

#include <boost/python.hpp>

#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
//#include <boost/python/args.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

namespace mmtbx { namespace rotamer {
namespace {

  boost::python::tuple
  getinitargs(xyzrad<> const& self)
  {
    return boost::python::make_tuple(self.sites_cart, self.radii);
  }

  void init_module()
  {
    using namespace boost::python;

    //
    typedef return_value_policy<return_by_value> rbv;
    class_<xyzrad<> >("xyzrad")
      .def(init<af::shared<scitbx::vec3<double> > const&,
                af::shared<double> const&,
                af::shared<double> const& >((arg("sites_cart"),
                                             arg("radii"),
                                             arg("weights"))))
      .add_property("sites_cart", make_getter(&xyzrad<>::sites_cart, rbv()))
      .add_property("radii",      make_getter(&xyzrad<>::radii, rbv()))
      .add_property("weights",    make_getter(&xyzrad<>::weights, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)

      .def(init<af::shared<scitbx::vec3<double> > const&,
                af::shared<double> const& >((arg("sites_cart"),
                                             arg("radii"))))
      .add_property("sites_cart", make_getter(&xyzrad<>::sites_cart, rbv()))
      .add_property("radii",    make_getter(&xyzrad<>::weights, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs)
    ;
    //

    class_<fit<> >("fit")

     .def(init<       double,
                      xyzrad<double> const&,
                      //af::shared<scitbx::vec3<double> >,
                      boost::python::list const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      af::const_ref<double, af::c_grid_padded<3> > const&,
                      //af::shared<scitbx::vec3<double> >,
                      xyzrad<double> const&,

                      cctbx::uctbx::unit_cell const&,
                      af::const_ref<std::size_t> const& ,
                      af::const_ref<std::size_t> const& ,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      double,
                      int >
                       ((arg("target_value"),
                         arg("xyzrad_bumpers"),
                         //arg("sites_cart_bumpers"),
                         arg("axes"),
                         arg("rotatable_points_indices"),
                         arg("angles_array"),
                         arg("density_map"),
                         arg("all_points"),
                         arg("unit_cell"),
                         arg("selection_clash"),
                         arg("selection_rsr"),
                         arg("sin_table"),
                         arg("cos_table"),
                         arg("step"),
                         arg("n"))))

     .def(init<double,
                      boost::python::list const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      af::const_ref<double, af::c_grid_padded<3> > const&,
                      af::shared<scitbx::vec3<double> >,
                      cctbx::uctbx::unit_cell const&,
                      af::const_ref<std::size_t> const& ,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      double,
                      int >
                       ((arg("target_value"),
                         arg("axes"),
                         arg("rotatable_points_indices"),
                         arg("angles_array"),
                         arg("density_map"),
                         arg("all_points"),
                         arg("unit_cell"),
                         arg("selection"),
                         arg("sin_table"),
                         arg("cos_table"),
                         arg("step"),
                         arg("n"))))

     .def(init<af::shared<scitbx::vec3<double> > const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      af::shared<scitbx::vec3<double> >,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      double,
                      int >
                       ((arg("sites_cart_start"),
                         arg("axes"),
                         arg("rotatable_points_indices"),
                         arg("angles_array"),
                         arg("all_points"),
                         arg("sin_table"),
                         arg("cos_table"),
                         arg("step"),
                         arg("n"))))

      .def("result", &fit<>::result)
      .def("score", &fit<>::score)
    ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::rotamer

BOOST_PYTHON_MODULE(mmtbx_rotamer_fit_ext)
{
  mmtbx::rotamer::init_module();
}
