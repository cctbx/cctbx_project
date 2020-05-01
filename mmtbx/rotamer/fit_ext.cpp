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
  getinitargs_moving(moving<> const& self)
  {
    return boost::python::make_tuple(
      self.sites_cart,
      self.sites_cart_start,
      self.radii,
      self.weights,
      self.bonded_pairs,
      self.max_map_value,
      self.min_map_value);
  }

  boost::python::tuple
  getinitargs_fixed(fixed<> const& self)
  {
    return boost::python::make_tuple(self.sites_cart, self.radii);
  }

  void init_module()
  {
    using namespace boost::python;

    //
    typedef return_value_policy<return_by_value> rbv;
    class_<moving<> >("moving")

      .def(init<af::shared<scitbx::vec3<double> > const&,
                af::shared<scitbx::vec3<double> > const&,
                af::shared<double> const&,
                af::shared<double> const&,
                boost::python::list const&,
                double const&,
                double const&>
                               ((arg("sites_cart"),
                                 arg("sites_cart_start"),
                                 arg("radii"),
                                 arg("weights"),
                                 arg("bonded_pairs"),
                                 arg("max_map_value"),
                                 arg("min_map_value"))))
      .add_property("sites_cart",       make_getter(&moving<>::sites_cart,       rbv()))
      .add_property("sites_cart_start", make_getter(&moving<>::sites_cart_start, rbv()))
      .add_property("radii",            make_getter(&moving<>::radii,            rbv()))
      .add_property("weights",          make_getter(&moving<>::weights,          rbv()))
      .add_property("bonded_pairs",     make_getter(&moving<>::bonded_pairs,     rbv()))
      .add_property("max_map_value",    make_getter(&moving<>::max_map_value,     rbv()))
      .add_property("min_map_value",    make_getter(&moving<>::min_map_value,     rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs_moving)
    ;

    //
    typedef return_value_policy<return_by_value> rbv;
    class_<fixed<> >("fixed")

      .def(init<af::shared<scitbx::vec3<double> > const&,
                af::shared<double> const& >((arg("sites_cart"),
                                             arg("radii"))))
      .add_property("sites_cart", make_getter(&fixed<>::sites_cart, rbv()))
      .add_property("radii",      make_getter(&fixed<>::radii, rbv()))
      .enable_pickling()
      .def("__getinitargs__", getinitargs_fixed)
    ;
    //

    class_<fit<> >("fit")

     .def(init<       fixed<double> const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      boost::python::list const&,
                      af::const_ref<double, af::c_grid_padded<3> > const&,
                      moving<double> const&,
                      cctbx::uctbx::unit_cell const&,
                      af::const_ref<std::size_t> const& ,
                      af::const_ref<std::size_t> const& ,
                      af::const_ref<double> const&,
                      af::const_ref<double> const&,
                      double,
                      int >
                       ((arg("fixed"),
                         arg("axes"),
                         arg("rotatable_points_indices"),
                         arg("angles_array"),
                         arg("density_map"),
                         arg("moving"),
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
      .def("score_start", &fit<>::score_start)
    ;

  }

} // namespace <anonymous>
}} // namespace mmtbx::rotamer

BOOST_PYTHON_MODULE(mmtbx_rotamer_fit_ext)
{
  mmtbx::rotamer::init_module();
}
