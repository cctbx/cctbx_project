#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/geometry_restraints/planarity.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct planarity_proxy_wrappers
  {
    typedef planarity_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("planarity_proxy", no_init)
        .def(init<
          af::shared<std::size_t> const&,
          af::shared<double> const&>(
            (arg_("i_seqs"), arg_("weights"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("weights", make_getter(&w_t::weights, rbv()))
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_planarity_proxy");
      }
    }
  };

  struct planarity_wrappers
  {
    typedef planarity w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_internal_reference<> rir;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("planarity", no_init)
        .def(init<
          af::shared<scitbx::vec3<double> > const&,
          af::shared<double> const&>(
            (arg_("sites"), arg_("weights"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  planarity_proxy const&>(
          (arg_("sites_cart"), arg_("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .add_property("weights", make_getter(&w_t::weights, rbv()))
        .def("deltas", &w_t::deltas, ccr())
        .def("rms_deltas", &w_t::rms_deltas)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
        .def("normal", &w_t::normal)
        .def("lambda_min", &w_t::lambda_min)
        .def("center_of_mass", &w_t::center_of_mass, ccr())
        .def("residual_tensor", &w_t::residual_tensor, ccr())
        .def("eigensystem", &w_t::eigensystem, rir())
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    planarity_proxy_wrappers::wrap();
    planarity_wrappers::wrap();
    def("planarity_deltas_rms", planarity_deltas_rms,
      (arg_("sites_cart"), arg_("proxies")));
    def("planarity_residuals", planarity_residuals,
      (arg_("sites_cart"), arg_("proxies")));
    def("planarity_residual_sum", planarity_residual_sum,
      (arg_("sites_cart"), arg_("proxies"), arg_("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_planarity() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
