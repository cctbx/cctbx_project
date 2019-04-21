#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <mmtbx/pair_interaction/interaction.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <boost/python.hpp>


namespace mmtbx { namespace pair_interaction {
  namespace bp = boost::python;

namespace {

  boost::python::tuple
  getinitargs(wfc<> const& self)
  {
    return boost::python::make_tuple(self.node_offsets,
      self.coefficients_of_first_derivative, self.coefficients_of_second_derivative);
  }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

   class_<wfc<> >("wfc")
     .def(init<af::shared<vec3<int> > const&,
               af::shared<vec3<int> > const&,
               af::shared<vec3<int> > const&,
               double const&,
               double const&,
               double const&,
               af::shared<vec3<double> > const&,
               af::shared<double> const&,
               int const&,
               double const& >(
                 (arg("node_offsets"),
                  arg("coefficients_of_first_derivative"),
                  arg("coefficients_of_second_derivative"),
                  arg("prefactor_of_first_derivative"),
                  arg("prefactor_of_second_derivative"),
                  arg("core_cutdens"),
                  arg("rr_array"),
                  arg("r_array"),
                  arg("ngrid"),
                  arg("zz")
                  )))
     .add_property("node_offsets", make_getter(&wfc<>::node_offsets, rbv()), make_setter(&wfc<>::node_offsets, dcp()))
     .add_property("coefficients_of_first_derivative", make_getter(&wfc<>::coefficients_of_first_derivative, rbv()), make_setter(&wfc<>::coefficients_of_first_derivative, dcp()))
     .add_property("coefficients_of_second_derivative", make_getter(&wfc<>::coefficients_of_second_derivative, rbv()), make_setter(&wfc<>::coefficients_of_second_derivative, dcp()))
     .add_property("prefactor_of_first_derivative", make_getter(&wfc<>::prefactor_of_first_derivative, rbv()), make_setter(&wfc<>::prefactor_of_first_derivative, dcp()))
     .add_property("prefactor_of_second_derivative", make_getter(&wfc<>::prefactor_of_second_derivative, rbv()), make_setter(&wfc<>::prefactor_of_second_derivative, dcp()))
     .add_property("a", make_getter(&wfc<>::a, rbv()), make_setter(&wfc<>::a, dcp()))
     .add_property("b", make_getter(&wfc<>::b, rbv()), make_setter(&wfc<>::b, dcp()))
     .add_property("position_max", make_getter(&wfc<>::position_max, rbv()), make_setter(&wfc<>::position_max, dcp()))
     .add_property("square_position_max", make_getter(&wfc<>::square_position_max, rbv()), make_setter(&wfc<>::square_position_max, dcp()))
     .add_property("ngrid", make_getter(&wfc<>::ngrid, rbv()), make_setter(&wfc<>::ngrid, dcp()))
     .add_property("grid_positions", make_getter(&wfc<>::grid_positions, rbv()), make_setter(&wfc<>::grid_positions, dcp()))
     //.add_property("grid_values", make_getter(&wfc<>::grid_values, rbv()), make_setter(&wfc<>::grid_values, dcp()))
     //.add_property("first_derivative_of_grid_values", make_getter(&wfc<>::first_derivative_of_grid_values, rbv()), make_setter(&wfc<>::first_derivative_of_grid_values, dcp()))
     //.add_property("second_derivative_of_grid_values", make_getter(&wfc<>::second_derivative_of_grid_values, rbv()), make_setter(&wfc<>::second_derivative_of_grid_values, dcp()))
     .add_property("core_cutdens", make_getter(&wfc<>::core_cutdens, rbv()), make_setter(&wfc<>::core_cutdens, dcp()))
     .enable_pickling()
     .def("__getinitargs__", getinitargs)
   ;

   class_<density_props<> >("density_props")
     .def(init<double,
               vec3<double> const&,
               mat3<double> const&
               >(
                 (arg("density"),
                  arg("gradient_vector"),
                  arg("hessian")
                  )))
     .def("add", &density_props<>::add, arg("density_props"))
     .def("has_silva_interaction", &density_props<>::has_silva_interaction)
     .add_property("density", make_getter(&density_props<>::density, rbv()), make_setter(&density_props<>::density, dcp()))
     .add_property("gradient_vector", make_getter(&density_props<>::gradient_vector, rbv()), make_setter(&density_props<>::gradient_vector, dcp()))
     .add_property("hessian", make_getter(&density_props<>::hessian, rbv()), make_setter(&density_props<>::hessian, dcp()))
     .add_property("gradient", make_getter(&density_props<>::gradient, rbv()), make_setter(&density_props<>::gradient, dcp()))
   ;

   def("hessian", (mat3<double>(*)(
                     vec3<double> const&,
                     double,
                     vec3<double> const&,
                     double,
                     double)) hessian, (arg("distanceVector"),
       arg("distanceReciprocal"),
       arg("distanceUnitVector"),arg("fac1"),arg("fac2")))
   ;

   def("atom_density_props", (density_props<double>(*)(
                     vec3<double> const&,
                     vec3<double> const&,
                     wfc<double>  const&)) atom_density_props, (
       arg("p"),
       arg("a_xyz"),
       arg("wfc_obj")))
   ;

   def("has_interaction_at_point", (bool(*)(
                     vec3<double> const&,
                     af::shared<vec3<double> > const&,
                     af::shared<int> const&,
                     //af::shared<wfc<double> >  const&
                     boost::python::list const&
                     )) has_interaction_at_point, (
       arg("p"),
       arg("a_xyz"),
       arg("element_flags"),
       arg("wfc_obj")
       ))
   ;

   def("points_and_pairs", (af::shared<vec3<int> >(*)(
                     vec3<int> const&,
                     double const&,
                     af::shared<vec3<double> > const&,
                     vec3<double> const&,
                     af::shared<int> const&,
                     af::shared<int> const&,
                     boost::python::list const&)) points_and_pairs, (
       arg("ngrid"),
       arg("step_size"),
       arg("xyz"),
       arg("xyz_min"),
       arg("atom_in_residue"),
       arg("element_flags"),
       arg("wfc_obj")))
   ;


  }

} // namespace <anonymous>
}} // namespace mmtbx::tls

BOOST_PYTHON_MODULE(mmtbx_pair_interaction_ext)
{
  mmtbx::pair_interaction::init_module();
}
