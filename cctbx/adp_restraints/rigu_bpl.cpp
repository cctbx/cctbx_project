#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/adp_restraints/rigu.h>

namespace cctbx { namespace adp_restraints {
namespace {

  struct rigu_proxy_wrappers
  {
    typedef rigu_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("rigu_proxy", no_init)
        .def(init<af::tiny<unsigned, 2> const&, double>((
           arg("i_seqs"),
           arg("weight"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_rigu_proxy")
        ;
      }
    }
  };

  struct rigu_wrappers
  {
    typedef rigu w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("rigu", no_init)
        .def(init<
           af::tiny<scitbx::vec3<double>, 2> const&,
           af::tiny<scitbx::sym_mat3<double>, 2> const&,
           double>(
          (arg("sites"),
           arg("u_cart"),
           arg("weight"))))
        .def(init<
          adp_restraint_params<double> const &,
           rigu_proxy const&>(
          (arg("params"),
           arg("proxy"))))
        .add_property("weight", make_getter(&w_t::weight, rbv()))
        .def("delta_33", &w_t::delta_33)
        .def("delta_13", &w_t::delta_13)
        .def("delta_23", &w_t::delta_23)
        .def("residual33", &w_t::residual33)
        .def("residual13", &w_t::residual13)
        .def("residual23", &w_t::residual23)
        .def("gradients33", &w_t::gradients33)
        .def("gradients13", &w_t::gradients13)
        .def("gradients23", &w_t::gradients23)
      ;
    }
  };

  void wrap_all() {
    using namespace boost::python;
    rigu_wrappers::wrap();

    rigu_proxy_wrappers::wrap();
      def("rigu_residual_sum",
        adp_restraint_residual_sum_aniso<rigu_proxy,rigu>::impl,
        (arg("params"),
         arg("proxies"),
         arg("gradients_aniso_cart")));
      def("rigu_residuals",
        adp_restraint_residuals<rigu_proxy,rigu>::impl,
        (arg("params"),
         arg("proxies")));
      def("rigu_deltas",
        rigu_deltas,
        (arg("params"),
         arg("proxies")));
  }

}

namespace boost_python {

  void wrap_rigu() { wrap_all(); }

}}}
