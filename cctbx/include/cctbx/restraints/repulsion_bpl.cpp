#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/restraints/repulsion.h>

namespace cctbx { namespace restraints {
namespace {

  struct repulsion_sym_proxy_wrappers
  {
    typedef repulsion_sym_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("repulsion_sym_proxy", no_init)
        .def(init<
          cctbx::crystal::direct_space_asu::asu_mapping_index_pair const&,
          double>(
            (arg_("pair"), arg_("vdw_radius"))))
        .def_readonly("pair", &w_t::pair)
        .def_readwrite("vdw_radius", &w_t::vdw_radius)
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<w_t, rir>::wrap(
          "shared_repulsion_sym_proxy");
      }
    }
  };

  struct repulsion_function_wrappers
  {
    typedef repulsion_function w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      class_<w_t>("repulsion_function", no_init)
        .def(init<optional<double, double, double, double> >(
          (arg_("c_rep"), arg_("k_rep"), arg_("irexp"), arg_("rexp"))))
        .def_readonly("c_rep", &w_t::c_rep)
        .def_readonly("k_rep", &w_t::k_rep)
        .def_readonly("irexp", &w_t::irexp)
        .def_readonly("rexp", &w_t::rexp)
      ;
    }
  };

  struct repulsion_wrappers
  {
    typedef repulsion w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("repulsion", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 2> const&,
                  double,
                  optional<repulsion_function const& > >(
          (arg_("sites"), arg_("vdw_radius"), arg_("function"))))
        .def_readonly("function", &w_t::function)
        .add_property("diff_vec", make_getter(&w_t::diff_vec, rbv()))
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  BOOST_PYTHON_FUNCTION_OVERLOADS(
    repulsion_deltas_overloads, repulsion_deltas, 3, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    repulsion_residuals_overloads, repulsion_residuals, 3, 4)
  BOOST_PYTHON_FUNCTION_OVERLOADS(
    repulsion_residual_sum_overloads, repulsion_residual_sum, 4, 6)

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    repulsion_sym_proxy_wrappers::wrap();
    repulsion_function_wrappers::wrap();
    repulsion_wrappers::wrap();
    def("repulsion_deltas", repulsion_deltas,
      repulsion_deltas_overloads(
        (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxies"),
         arg_("function"))));
    def("repulsion_residuals", repulsion_residuals,
      repulsion_residuals_overloads(
        (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxies"),
         arg_("function"))));
    def("repulsion_residual_sum", repulsion_residual_sum,
      repulsion_residual_sum_overloads(
        (arg_("sites_cart"), arg_("asu_mappings"), arg_("proxies"),
         arg_("gradient_array"), arg_("function"), arg_("disable_cache"))));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_repulsion() { wrap_all(); }

}}} // namespace cctbx::boost_python
