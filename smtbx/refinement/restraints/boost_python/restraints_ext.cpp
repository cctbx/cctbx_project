#include <boost/python/module.hpp>
#include <boost/python/class.hpp>
#include <boost/python/def.hpp>

#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>

#include <cctbx/geometry_restraints/angle.h>
#include <cctbx/geometry_restraints/bond.h>
#include <cctbx/geometry_restraints/dihedral.h>
#include <cctbx/geometry_restraints/bond_similarity.h>
#include <cctbx/adp_restraints/adp_similarity.h>
#include <cctbx/adp_restraints/rigid_bond.h>
#include <cctbx/adp_restraints/isotropic_adp.h>
#include <cctbx/adp_restraints/fixed_u_eq_adp.h>
#include <cctbx/restraints.h>

#include <smtbx/import_scitbx_af.h>

namespace smtbx { namespace refinement { namespace restraints {

namespace boost_python {

  struct functions_wrapper {
    template <typename ProxyType, typename RestraintType>
    static void wrap() {
      using namespace boost::python;
      def("linearise_restraints",
      (void(*) (
        cctbx::uctbx::unit_cell const &,
        cctbx::adp_restraints::adp_restraint_params<double> const &,
        cctbx::xray::parameter_map<cctbx::xray::scatterer<double> > const &,
        af::const_ref<ProxyType> const &,
        cctbx::restraints::linearised_eqns_of_restraint<double> &))
        cctbx::restraints::linearise_restraints<
          double, cctbx::adp_restraints::adp_restraint_params<double>,
          ProxyType, RestraintType>,
        (arg("unit_cell"),
         arg("params"),
         arg("parameter_map"),
         arg("proxies"),
         arg("linearised_eqns")));
    }
  };

  template <typename FloatType>
  struct linearised_eqns_of_restraint_wrapper
  {
    typedef cctbx::restraints::linearised_eqns_of_restraint<FloatType> wt;

    static void wrap(char const *name) {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<wt>(name, no_init)
        .def(init<std::size_t, std::size_t>
             ((arg("n_restraints"),
               arg("n_crystallographic_params"))))
        .def_readwrite("design_matrix", &wt::design_matrix)
        .add_property("deltas", make_getter(&wt::deltas, rbv()))
        .add_property("weights", make_getter(&wt::weights, rbv()))
        .def("n_crystallographic_params", &wt::n_crystallographic_params)
        .def("n_restraints", &wt::n_restraints)
        .def("add_equation", &wt::add_equation)
        ;
    }
  };

  template <typename FloatType, typename ProxyType, typename RestraintType>
  struct geom_res_linearise_restraints_wrapper
  {
    static void wrap() {
      using namespace boost::python;
      def("linearise_restraints",
        (void(*) (
          cctbx::uctbx::unit_cell const &,
          af::const_ref<scitbx::vec3<FloatType> > const &,
          cctbx::xray::parameter_map<cctbx::xray::scatterer<FloatType> > const &,
          af::const_ref<ProxyType> const &,
          cctbx::restraints::linearised_eqns_of_restraint<FloatType> &))
          cctbx::restraints::linearise_restraints<
            FloatType, ProxyType, RestraintType>,
          (arg("unit_cell"),
           arg("sites_cart"),
           arg("parameter_map"),
           arg("proxies"),
           arg("restraints_matrix")));
    }
  };

  namespace geom_res = cctbx::geometry_restraints;
  namespace adp_res = cctbx::adp_restraints;

  void wrap_least_squares_restraints() {
    using namespace boost::python;

    linearised_eqns_of_restraint_wrapper<
      double>::wrap("linearised_eqns_of_restraint");

    geom_res_linearise_restraints_wrapper<
      double, geom_res::angle_proxy, geom_res::angle>::wrap();
    geom_res_linearise_restraints_wrapper<
      double, geom_res::bond_simple_proxy, geom_res::bond>::wrap();
    geom_res_linearise_restraints_wrapper<
      double, geom_res::dihedral_proxy, geom_res::dihedral>::wrap();
    //geom_res_linearise_restraints_wrapper<
    //  double, double, geom_res::planarity_proxy, geom_res::planarity>::wrap();
    geom_res_linearise_restraints_wrapper<
      double, geom_res::bond_similarity_proxy, geom_res::bond_similarity>::wrap();
    // wrap linearisation functions
    functions_wrapper::wrap<
      adp_res::isotropic_adp_proxy, adp_res::isotropic_adp>();
    functions_wrapper::wrap<
      adp_res::fixed_u_eq_adp_proxy, adp_res::fixed_u_eq_adp>();
    functions_wrapper::wrap<
      adp_res::adp_similarity_proxy, adp_res::adp_similarity>();
    functions_wrapper::wrap<
      adp_res::adp_u_eq_similarity_proxy, adp_res::adp_u_eq_similarity>();
    functions_wrapper::wrap<
      adp_res::rigid_bond_proxy, adp_res::rigid_bond>();
    functions_wrapper::wrap<
      adp_res::adp_volume_similarity_proxy, adp_res::adp_volume_similarity>();
  }

  void wrap_origin_fixing_restraints();

  namespace {
    void init_module() {
      wrap_least_squares_restraints();
      wrap_origin_fixing_restraints();
    }
  }

}}}} // namespace smtbx::refinement::restraints::boost_python

BOOST_PYTHON_MODULE(smtbx_refinement_restraints_ext)
{
  smtbx::refinement::restraints::boost_python::init_module();
}
