#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/fixed_u_eq_adp.h>
#include <scitbx/boost_python/container_conversions.h>


namespace cctbx { namespace adp_restraints {

namespace {

  struct fixed_u_eq_adp_proxy_wrappers  {
    typedef fixed_u_eq_adp_proxy w_t;

    static void wrap() {
      using namespace boost::python;
      class_<w_t, bases<adp_restraint_proxy<1> > >
        ("fixed_u_eq_adp_proxy", no_init)
        .def(init<
              af::tiny<unsigned, 1> const &,
              double,
              double>(
          (arg("i_seqs"),
           arg("weight"),
           arg("u_eq_ideal"))))
        .def_readonly("u_eq_ideal", &w_t::u_eq_ideal)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_fixed_u_eq_adp_proxy");
      }
    }
  };

  struct fixed_u_eq_adp_wrappers  {
    typedef fixed_u_eq_adp w_t;

    static void wrap() {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;

      class_<w_t, bases<adp_restraint_base_1<1> > >
            ("fixed_u_eq_adp", no_init)
        .def(init<
            scitbx::sym_mat3<double> const &,
            double,
            double>(
          (arg("u_cart"),
           arg("weight"),
           arg("u_eq_ideal"))))
        .def(init<
            double,
            double,
            double>(
          (arg("u_iso"),
           arg("weight"),
           arg("u_eq_ideal"))))
        .def(init<
            adp_restraint_params<double> const &,
            fixed_u_eq_adp_proxy const &>(
          (arg("params"),
           arg("proxy"))))
        .def_readonly("u_eq_ideal", &w_t::u_eq_ideal)
      ;
    }
  };

  void wrap_all() {
    using namespace boost::python;
    fixed_u_eq_adp_wrappers::wrap();
    fixed_u_eq_adp_proxy_wrappers::wrap();
  }

}

namespace boost_python {

  void wrap_fixed_u_eq_adp() { wrap_all(); }

}}}
