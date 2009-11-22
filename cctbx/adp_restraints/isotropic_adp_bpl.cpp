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
      class_<w_t>("isotropic_adp_proxy", no_init)
        .def(init<unsigned, double>((
           arg("i_seq"),
           arg("weight"))))
        .add_property("i_seq", make_getter(&w_t::i_seq, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
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
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("isotropic_adp", no_init)
        .def(init<
           scitbx::sym_mat3<double> const&,
           double>(
          (arg("u_cart"),
           arg("weight"))))
        .def(init<
           af::const_ref<scitbx::sym_mat3<double> > const&,
           isotropic_adp_proxy const&>(
          (arg("u_cart"),
           arg("proxy"))))
        .add_property("u_cart", make_getter(&w_t::u_cart, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
        .def("deltas", &w_t::deltas)
        .def("rms_deltas", &w_t::rms_deltas)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    isotropic_adp_wrappers::wrap();
    isotropic_adp_proxy_wrappers::wrap();
    def("isotropic_adp_residual_sum", isotropic_adp_residual_sum,
      (arg("u_cart"),
       arg("proxies"),
       arg("gradients_aniso_cart")));
    def("isotropic_adp_residuals", isotropic_adp_residuals,
      (arg("u_cart"),
       arg("proxies")));
    def("isotropic_adp_deltas_rms", isotropic_adp_deltas_rms,
      (arg("u_cart"),
       arg("proxies")));
  }

}

namespace boost_python {

  void
  wrap_isotropic_adp() { wrap_all(); }

}}}
