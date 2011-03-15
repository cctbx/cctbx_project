
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <cctbx/adp_restraints/aniso_restraints.h>

namespace cctbx {namespace adp_restraints {
namespace {

  void wrap_all ()
  {
    using namespace boost::python;
    using namespace cctbx::geometry_restraints;
    using namespace cctbx::xray;
    typedef eval_adp_aniso_restraints w_t;
    typedef return_value_policy<return_by_value> rbv;
    class_<w_t>("eval_adp_aniso_restraints", no_init)
      .def(init<
        af::const_ref<scatterer<> > const&,
        af::const_ref<scitbx::sym_mat3<double> > const&,
        af::const_ref<double> const&,
        af::const_ref<bond_simple_proxy> const&,
        af::const_ref<bool> const&,
        af::const_ref<bool> const&,
        unsigned,
        bool>((
          arg("scatterers"),
          arg("u_cart"),
          arg("u_iso"),
          arg("bond_proxies"),
          arg("selection"),
          arg("hd_selection"),
          arg("n_grad_u_iso"),
          arg("use_hd"))))
      .def("gradients_iso", &w_t::gradients_iso)
      .def("gradients_aniso_cart", &w_t::gradients_aniso_cart)
      .def_readonly("target", &w_t::target)
      .def_readonly("number_of_restraints", &w_t::number_of_restraints)
    ;
  }
}
namespace boost_python {
  void wrap_aniso_restraints () { wrap_all(); }
}
}}
