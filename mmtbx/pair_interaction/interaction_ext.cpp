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
  getinitargs(wfc const& self)
  {
    return boost::python::make_tuple(self.ngrid, self.zz, self.r_array, self.wfcin_array, self.occ_electrons);
  }

  void init_module()
  {
    using namespace boost::python;
    using boost::python::arg;
    typedef return_value_policy<return_by_value> rbv;
    typedef default_call_policies dcp;

   class_<wfc>("wfc")
     .def(init<int const&,
               double const&,
               af::shared<double>  const&,
               boost::python::list const&,
               af::shared<double>
               >(
                 (arg("ngrid"),
                  arg("zz"),
                  arg("r_array"),
                  arg("wfcin_array"),
                  arg("occ_electrons")
                  )))
     .add_property("ngrid", make_getter(&wfc::ngrid, rbv()), make_setter(&wfc::ngrid, dcp()))
     .add_property("zz", make_getter(&wfc::zz, rbv()), make_setter(&wfc::zz, dcp()))
     .add_property("r_array", make_getter(&wfc::r_array, rbv()), make_setter(&wfc::r_array, dcp()))
     .add_property("wfcin_array", make_getter(&wfc::wfcin_array, rbv()), make_setter(&wfc::wfcin_array, dcp()))
     .add_property("occ_electrons", make_getter(&wfc::occ_electrons, rbv()), make_setter(&wfc::occ_electrons, dcp()))

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
     .def("get_dori_value", &density_props<>::get_dori_value)
     .def("get_sedd_value", &density_props<>::get_sedd_value)
     .def("cal_silva", &density_props<>::cal_silva)
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
                     vec3<double> const& ,
                     wfc          const&
                     )) atom_density_props, (
       arg("p"),
       arg("a_xyz"),
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
                     boost::python::list const&,
                     std::string const & silva_type)) points_and_pairs, (
       arg("ngrid"),
       arg("step_size"),
       arg("xyz"),
       arg("xyz_min"),
       arg("atom_in_residue"),
       arg("element_flags"),
       arg("wfc_obj"),
       arg("silva_type")
       ))
   ;


  }

} // namespace <anonymous>
}} // namespace mmtbx::tls

BOOST_PYTHON_MODULE(mmtbx_pair_interaction_ext)
{
  mmtbx::pair_interaction::init_module();
}
