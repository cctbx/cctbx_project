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
#include <cctbx/geometry_restraints/chirality.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct chirality_proxy_wrappers
  {
    typedef chirality_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("chirality_proxy", no_init)
        .def(init<af::tiny<unsigned, 4> const&, double, bool, double>(
          (arg_("i_seqs"),
           arg_("volume_ideal"),
           arg_("both_signs"),
           arg_("weight"))))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readonly("volume_ideal", &w_t::volume_ideal)
        .def_readonly("both_signs", &w_t::both_signs)
        .def_readonly("weight", &w_t::weight)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_chirality_proxy");
      }
    }
  };

  struct chirality_wrappers
  {
    typedef chirality w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("chirality", no_init)
        .def(init<af::tiny<scitbx::vec3<double>, 4> const&,
                  double, bool, double>(
          (arg_("sites"),
           arg_("volume_ideal"),
           arg_("both_signs"),
           arg_("weight"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  chirality_proxy const&>(
          (arg_("sites_cart"), arg_("proxy"))))
        .add_property("sites", make_getter(&w_t::sites, rbv()))
        .def_readonly("volume_ideal", &w_t::volume_ideal)
        .def_readonly("both_signs", &w_t::both_signs)
        .def_readonly("weight", &w_t::weight)
        .def_readonly("volume_model", &w_t::volume_model)
        .def_readonly("delta_sign", &w_t::delta_sign)
        .def_readonly("delta", &w_t::delta)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    chirality_proxy_wrappers::wrap();
    chirality_wrappers::wrap();
    def("chirality_deltas", chirality_deltas,
      (arg_("sites_cart"), arg_("proxies")));
    def("chirality_residuals", chirality_residuals,
      (arg_("sites_cart"), arg_("proxies")));
    def("chirality_residual_sum", chirality_residual_sum,
      (arg_("sites_cart"), arg_("proxies"), arg_("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_chirality() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
