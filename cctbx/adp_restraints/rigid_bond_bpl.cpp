#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/rigid_bond.h>

namespace cctbx { namespace adp_restraints {
namespace {

  struct rigid_bond_pair_wrappers
  {

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_;
      class_<rigid_bond_pair>("rigid_bond_pair",
                               init<vec3<double> const&,
                                    vec3<double> const&,
                                    sym_mat3<double> const&,
                                    sym_mat3<double> const&,
                                    cctbx::uctbx::unit_cell const&>())
        .def("z_12", &rigid_bond_pair::z_12)
        .def("z_21", &rigid_bond_pair::z_21)
        .def("delta_z", &rigid_bond_pair::delta_z)
      ;
    }
  };

  struct rigid_bond_proxy_wrappers
  {
    typedef rigid_bond_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("rigid_bond_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg("i_seqs"),
           arg("weight"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_rigid_bond_proxy")
        ;
      }
    }
  };

  struct rigid_bond_wrappers
  {
    typedef rigid_bond w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("rigid_bond", no_init)
        .def(init<
           af::tiny<scitbx::vec3<double>, 2> const&,
           af::tiny<scitbx::sym_mat3<double>, 2> const&,
           double>(
          (arg("sites"),
           arg("u_cart"),
           arg("weight"))))
        .def(init<
          adp_restraint_params<double> const &,
           rigid_bond_proxy const&>(
          (arg("params"),
           arg("proxy"))))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
        .def("z_12", &w_t::z_12)
        .def("z_21", &w_t::z_21)
        .def("delta_z", &w_t::delta_z)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  void wrap_all() {
    using namespace boost::python;
    rigid_bond_pair_wrappers::wrap();
    rigid_bond_wrappers::wrap();

    rigid_bond_proxy_wrappers::wrap();
      def("rigid_bond_residual_sum",
        adp_restraint_residual_sum_aniso<rigid_bond_proxy,rigid_bond>::impl,
        (arg("params"),
         arg("proxies"),
         arg("gradients_aniso_cart")));
      def("rigid_bond_residuals",
        adp_restraint_residuals<rigid_bond_proxy,rigid_bond>::impl,
        (arg("params"),
         arg("proxies")));
      def("rigid_bond_deltas",
        rigid_bond_deltas,
        (arg("params"),
         arg("proxies")));
  }

}

namespace boost_python {

  void wrap_rigid_bond() { wrap_all(); }

}}}
