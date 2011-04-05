#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>
#include <cctbx/geometry_restraints/nonbonded_sorted.h>

namespace cctbx { namespace geometry_restraints {
namespace {

  struct nonbonded_sorted_asu_proxies_base_wrappers
  {
    typedef nonbonded_sorted_asu_proxies_base w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      typedef return_value_policy<copy_const_reference> ccr;
      class_<w_t>("nonbonded_sorted_asu_proxies_base", no_init)
        .def(init<
          boost::shared_ptr<asu_mappings> const&>(
            (arg("asu_mappings"))))
        .def("asu_mappings", &w_t::asu_mappings, ccr())
        .def("process",
          (void(w_t::*)(nonbonded_simple_proxy const&)) &w_t::process,
            (arg("proxy")))
        .def("process",
          (void(w_t::*)(af::const_ref<nonbonded_simple_proxy> const&))
            &w_t::process,
          (arg("proxies")))
        .def("process",
          (void(w_t::*)(nonbonded_asu_proxy const&, bool)) &w_t::process,
            (arg("proxy"), arg("sym_excl_flag")=false))
        .def("process",
          (void(w_t::*)(af::const_ref<nonbonded_asu_proxy> const&))
            &w_t::process,
          (arg("proxies")))
        .def("n_total", &w_t::n_total)
        .def_readonly("simple", &w_t::simple)
        .def_readonly("asu", &w_t::asu)
      ;
    }
  };

  struct nonbonded_sorted_asu_proxies_wrappers
  {
    typedef nonbonded_sorted_asu_proxies w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t, bases<nonbonded_sorted_asu_proxies_base> >(
        "nonbonded_sorted_asu_proxies", no_init)
        .def(init<
          boost::shared_ptr<asu_mappings> const&>(
            (arg("asu_mappings"))))
        .def(init<
          af::const_ref<std::size_t> const&,
          af::const_ref<std::size_t> const&,
          af::const_ref<std::size_t> const&,
          af::const_ref<std::size_t> const&,
          nonbonded_params const&,
          af::const_ref<std::string> const&,
          double,
          double,
          std::vector<crystal::pair_asu_table<> > const&>((
            arg("model_indices"),
            arg("conformer_indices"),
            arg("sym_excl_indices"),
            arg("donor_acceptor_excl_groups"),
            arg("nonbonded_params"),
            arg("nonbonded_types"),
            arg("nonbonded_distance_cutoff_plus_buffer"),
            arg("min_cubicle_edge"),
            arg("shell_asu_tables"))))
        .def_readonly("n_unknown_nonbonded_type_pairs",
          &w_t::n_unknown_nonbonded_type_pairs)
        .def_readonly("min_vdw_distance", &w_t::min_vdw_distance)
        .def_readonly("max_vdw_distance", &w_t::max_vdw_distance)
      ;
    }
  };

  void
  wrap_all()
  {
    nonbonded_sorted_asu_proxies_base_wrappers::wrap();
    nonbonded_sorted_asu_proxies_wrappers::wrap();
  }

} // namespace <anonymous>

namespace boost_python {

  void
  wrap_nonbonded_sorted() { wrap_all(); }

}}} // namespace cctbx::geometry_restraints::boost_python
