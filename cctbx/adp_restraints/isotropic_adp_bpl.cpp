#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/isotropic_adp.h>
#include <scitbx/boost_python/container_conversions.h>


namespace cctbx { namespace adp_restraints {

namespace {


  struct isotropic_adp_proxy_wrappers
  {
    typedef isotropic_adp_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<adp_restraint_proxy<1> > >
        ("isotropic_adp_proxy", no_init)
        .def(init<
           af::tiny<unsigned, 1> const &,
           double>(
          (arg("i_seqs"),
           arg("weight"))))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_isotropic_adp_proxy")
        ;
      }
    }
  };

  struct isotropic_adp_wrappers
  {
    typedef isotropic_adp w_t;

    static void
    wrap() {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t, bases<adp_restraint_base_6<1> > >
            ("isotropic_adp", no_init)
        .def(init<
            scitbx::sym_mat3<double> const &,
            double>(
          (arg("u_cart"),
           arg("weight"))))
        .def(init<
            adp_restraint_params<double> const &,
            isotropic_adp_proxy const &>(
          (arg("params"),
           arg("proxy"))))
      ;
    }
  };

  void wrap_all() {
    using namespace boost::python;
    isotropic_adp_wrappers::wrap();
    isotropic_adp_proxy_wrappers::wrap();
  }

}

namespace boost_python {

  void
  wrap_isotropic_adp() { wrap_all(); }

}}}
