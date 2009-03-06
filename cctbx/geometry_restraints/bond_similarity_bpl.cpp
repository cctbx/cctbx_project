#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <boost/python/return_by_value.hpp>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <scitbx/array_family/selections.h>
#include <scitbx/stl/map_wrapper.h>
#include <cctbx/geometry_restraints/bond_similarity.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct bond_similarity_proxy_wrappers
  {
    typedef bond_similarity_proxy w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("bond_similarity_proxy", no_init)
        .def(init<
          af::shared<af::tiny<unsigned, 2> >,
          af::shared<sgtbx::rt_mx>,
          af::shared<double> const&>((
            arg_("i_seqs"),
            arg_("sym_ops"),
            arg_("weights"))))
        .add_property("i_seqs", make_getter(&w_t::i_seqs, rbv()))
        .add_property("weights", make_getter(&w_t::weights, rbv()))
        .add_property("sym_ops", make_getter(&w_t::sym_ops, rbv()))
      ;
      {
        typedef return_internal_reference<> rir;
        scitbx::af::boost_python::shared_wrapper<
          bond_similarity_proxy, rir>::wrap(
          "shared_bond_similarity_proxy");
      }
    }
  };

  struct bond_similarity_wrappers
  {
    typedef bond_similarity w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      typedef return_value_policy<return_by_value> rbv;
      scitbx::af::boost_python::shared_wrapper<
        af::tiny<scitbx::vec3<double>, 2>, rbv>::wrap(
        "sites_array");
      class_<w_t>("bond_similarity", no_init)
        .def(init<
          af::shared<af::tiny<scitbx::vec3<double>, 2> > const&,
          af::shared<double> const&>((
            arg_("sites_array"),
            arg_("weights"))))
        .def(init<uctbx::unit_cell const&,
                  af::const_ref<scitbx::vec3<double> > const&,
                  bond_similarity_proxy const&>(
          (arg_("unit_cell"), arg_("sites_cart"), arg_("proxy"))))
        .add_property("sites_array", make_getter(
          &w_t::sites_array, rbv()))
        .add_property("weights", make_getter(&w_t::weights, rbv()))
        .def("deltas", &w_t::deltas, ccr())
        .def("rms_deltas", &w_t::rms_deltas)
        .def("residual", &w_t::residual)
        .def("gradients", &w_t::gradients)
        .def("mean_distance", &w_t::mean_distance)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    bond_similarity_proxy_wrappers::wrap();
    bond_similarity_wrappers::wrap();
    def("bond_similarity_deltas_rms",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_similarity_proxy> const&))
      bond_similarity_deltas_rms,
      (arg_("unit_cell"), arg_("sites_cart"), arg_("proxies")));
    def("bond_similarity_residuals",
      (af::shared<double>(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_similarity_proxy> const&))
      bond_similarity_residuals,
      (arg_("unit_cell"), arg_("sites_cart"), arg_("proxies")));
   def("bond_similarity_residual_sum",
      (double(*)(
        uctbx::unit_cell const&,
        af::const_ref<scitbx::vec3<double> > const&,
        af::const_ref<bond_similarity_proxy> const&,
        af::ref<scitbx::vec3<double> > const&))
      bond_similarity_residual_sum,
      (arg_("unit_cell"),
       arg_("sites_cart"),
       arg_("proxies"),
       arg_("gradient_array")));
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_bond_similarity() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
