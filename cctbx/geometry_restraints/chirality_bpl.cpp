#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_internal_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/boost_python/is_polymorphic_workaround.h>
#include <cctbx/geometry_restraints/chirality.h>
#include <cctbx/geometry_restraints/proxy_select.h>

SCITBX_BOOST_IS_POLYMORPHIC_WORKAROUND(cctbx::geometry_restraints::chirality)

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
        .def(init<af::tiny<unsigned, 4> const&, double, bool, double>((
          arg("i_seqs"),
          arg("volume_ideal"),
          arg("both_signs"),
          arg("weight"))))
        .def(init<af::tiny<unsigned, 4> const&, w_t const&>((
          arg("i_seqs"),
          arg("proxy"))))
        .def("scale_weight", &w_t::scale_weight, (arg("factor")))
        .def("sort_i_seqs", &w_t::sort_i_seqs)
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .def_readonly("volume_ideal", &w_t::volume_ideal)
        .def_readonly("both_signs", &w_t::both_signs)
        .def_readwrite("weight", &w_t::weight)
      ;
      {
        scitbx::af::boost_python::shared_wrapper<w_t>::wrap(
          "shared_chirality_proxy")
          .def("proxy_select",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              std::size_t,
              af::const_ref<std::size_t> const&))
                shared_proxy_select, (
            arg("n_seq"), arg("iselection")))
          .def("proxy_remove",
            (af::shared<w_t>(*)(
              af::const_ref<w_t> const&,
              af::const_ref<bool> const&))
                shared_proxy_remove, (
            arg("selection")))
        ;
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
          (arg("sites"),
           arg("volume_ideal"),
           arg("both_signs"),
           arg("weight"))))
        .def(init<af::const_ref<scitbx::vec3<double> > const&,
                  chirality_proxy const&>(
          (arg("sites_cart"), arg("proxy"))))
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
      (arg("sites_cart"), arg("proxies")));
    def("chirality_residuals", chirality_residuals,
      (arg("sites_cart"), arg("proxies")));
    def("chirality_residual_sum", chirality_residual_sum,
      (arg("sites_cart"), arg("proxies"), arg("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_chirality() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
