#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/restraints/pair_proxies.h>

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
      class_<w_t>("pair_proxies", no_init)
        .def(init<
          af::const_ref<bond_params_dict> const&,
          af::const_ref<std::string> const&,
          repulsion_distance_table const&,
          repulsion_radius_table const&,
          double,
          std::vector<crystal::pair_asu_table<> > const&,
          double,
          double,
          double,
          double>((
           arg_("bond_params_table"),
           arg_("repulsion_types"),
           arg_("repulsion_distance_table"),
           arg_("repulsion_radius_table"),
           arg_("repulsion_distance_default"),
           arg_("shell_asu_tables"),
           arg_("bonded_distance_cutoff"),
           arg_("nonbonded_distance_cutoff"),
           arg_("nonbonded_buffer"),
           arg_("vdw_1_4_factor"))))
        .def_readonly("bond_proxies", &w_t::bond_proxies)
        .def_readonly("repulsion_proxies", &w_t::repulsion_proxies)
        .def_readonly("n_bonded", &w_t::n_bonded)
        .def_readonly("n_1_3", &w_t::n_1_3)
        .def_readonly("n_1_4", &w_t::n_1_4)
        .def_readonly("n_nonbonded", &w_t::n_nonbonded)
        .def_readonly("n_unknown_repulsion_type_pairs",
          &w_t::n_unknown_repulsion_type_pairs)
      ;
    }
  };

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    def("add_pairs", add_pairs, (
      arg_("pair_asu_table"), arg_("bond_simple_proxies")));
    pair_proxies_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_pair_proxies() { wrap_all(); }

}}} // namespace cctbx::restraints::boost_python
