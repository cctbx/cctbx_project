#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <cctbx/restraints/pair_proxies.h>
#include <cctbx/xray/scatterer.h>

namespace cctbx { namespace restraints {
namespace {

  struct pair_proxies_wrappers
  {
    typedef pair_proxies w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef boost::python::arg arg_; // gcc 2.96 workaround
      typedef return_value_policy<return_by_value> rbv;
      class_<w_t>("pair_proxies", no_init)
        .def(init<
          af::const_ref<xray::scatterer<> > const&,
          restraints::bond_params_table const&,
          restraints::repulsion_distance_table const&,
          std::vector<crystal::pair_asu_table<> > const&,
          af::const_ref<double> const&,
          double,
          double,
          double>((
           arg_("scatterers"),
           arg_("bond_params_table"),
           arg_("repulsion_distance_table"),
           arg_("shell_asu_tables"),
           arg_("shell_distance_cutoffs"),
           arg_("nonbonded_distance_cutoff"),
           arg_("nonbonded_buffer"),
           arg_("vdw_1_4_factor"))))
        .add_property("bond_asu_proxies",
          make_getter(&w_t::bond_asu_proxies, rbv()))
        .add_property("repulsion_asu_proxies",
          make_getter(&w_t::repulsion_asu_proxies, rbv()))
      ;
    }
  };

  void
  wrap_all()
  {
    pair_proxies_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_pair_proxies() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
