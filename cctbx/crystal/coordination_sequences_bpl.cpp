#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/boost_python/container_conversions.h>
#include <cctbx/crystal/coordination_sequences.h>
#include <cctbx/crystal/workarounds_bpl.h>

namespace cctbx { namespace crystal { namespace coordination_sequences {
namespace {

  void
  wrap_all()
  {
    using namespace boost::python;
    def("coordination_sequences_simple", coordination_sequences::simple, (
      arg("pair_asu_table"),
      arg("max_shell")));
    def("coordination_sequences_simple_sym",
         coordination_sequences::simple_sym, (
      arg("full_pair_sym_table"),
      arg("site_symmetry_table"),
      arg("max_shell")));
    def("coordination_sequences_shell_asu_tables",
      coordination_sequences::shell_asu_tables, (
        arg("pair_asu_table"),
        arg("max_shell")));
    def("coordination_sequences_shell_sym_tables",
      coordination_sequences::shell_sym_tables, (
        arg("full_pair_sym_table"),
        arg("site_symmetry_table"),
        arg("max_shell")));
    {
      using namespace scitbx::boost_python::container_conversions;
      tuple_mapping<
        std::vector<pair_asu_table<> >, variable_capacity_policy>();
      tuple_mapping<
        std::vector<pair_sym_table>, variable_capacity_policy>();
    }
  }

}} // namespace coordination_sequences::<anonymous>

namespace boost_python {

  void
  wrap_coordination_sequences() { coordination_sequences::wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
