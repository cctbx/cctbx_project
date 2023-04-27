#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/other_restraints/sump.h>

namespace cctbx { namespace other_restraints {
namespace {

  struct sump_proxy_wrapper
  {
    typedef sump_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<w_t>("sump_proxy", no_init)
        .def(init<af::shared<unsigned> const&,
                  af::shared<double> const& , double, double>((
           arg("i_seqs"),
           arg("coefficients"),
           arg("weight"),
           arg("target"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv))
        .add_property("coefficients", make_getter(&w_t::coefficients, rbv))
        .add_property("weight", &w_t::weight)
        .add_property("target", &w_t::target)
        ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_sump_proxy")
        ;
      }
    }
  };

  struct sump_wrapper {
    typedef sump w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      return_value_policy<return_by_value> rbv;
      class_<w_t>("sump", no_init)
        .def(init<
           af::shared<double> const&,
          af::shared<double> const&,
          double, double>(
          (arg("occupancies"),
           arg("coefficients"),
           arg("weight"), arg("target"))))
        .def(init<
          af::shared<double> const&,
          sump_proxy const&>(
            (arg("occupancies"),
             arg("proxy"))))
        .add_property("weight", &w_t::weight)
        .add_property("target", &w_t::target)
        .add_property("delta", &w_t::delta)
        .add_property("residual", &w_t::residual)
        .add_property("coefficients", make_getter(&w_t::coefficients, rbv))
        .def("linearise", &w_t::linearise)
      ;
    }
  };

  void wrap_all() {
    using namespace boost::python;
    sump_proxy_wrapper::wrap();
    sump_wrapper::wrap();

      //def("rigid_bond_residual_sum",
      //  adp_restraint_residual_sum_aniso<rigid_bond_proxy,rigid_bond>::impl,
      //  (arg("params"),
      //   arg("proxies"),
      //   arg("gradients_aniso_cart")));
      //def("rigid_bond_residuals",
      //  adp_restraint_residuals<rigid_bond_proxy,rigid_bond>::impl,
      //  (arg("params"),
      //   arg("proxies")));
      //def("rigid_bond_deltas",
      //  rigid_bond_deltas,
      //  (arg("params"),
      //   arg("proxies")));
  }

}

namespace boost_python {

  void wrap_sump() { wrap_all(); }

}}}
