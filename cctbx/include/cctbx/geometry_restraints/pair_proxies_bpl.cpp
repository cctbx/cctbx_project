#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/geometry_restraints/pair_proxies.h>

namespace cctbx { namespace geometry_restraints {
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
          af::const_ref<bond_params_dict> const&>((
           arg_("bond_params_table"))))
        .def(init<
          nonbonded_params const&,
          af::const_ref<std::string> const&,
          af::const_ref<bond_params_dict> const&,
          std::vector<crystal::pair_asu_table<> > const&,
          double,
          double,
          double>((
           arg_("nonbonded_params"),
           arg_("nonbonded_types"),
           arg_("bond_params_table"),
           arg_("shell_asu_tables"),
           arg_("bonded_distance_cutoff"),
           arg_("nonbonded_distance_cutoff"),
           arg_("nonbonded_buffer"))))
        .def_readonly("bond_proxies", &w_t::bond_proxies)
        .def_readonly("nonbonded_proxies", &w_t::nonbonded_proxies)
        .def_readonly("n_bonded", &w_t::n_bonded)
        .def_readonly("n_1_3", &w_t::n_1_3)
        .def_readonly("n_1_4", &w_t::n_1_4)
        .def_readonly("n_nonbonded", &w_t::n_nonbonded)
        .def_readonly("n_unknown_nonbonded_type_pairs",
          &w_t::n_unknown_nonbonded_type_pairs)
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

}}} // namespace cctbx::geometry_restraints::boost_python
