#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <cctbx/crystal/coordination_sequences.h>

namespace cctbx { namespace crystal { namespace coordination_sequences {
namespace {

  void
  wrap_all()
  {
    using namespace boost::python;
    typedef boost::python::arg arg_; // gcc 2.96 workaround
    def("coordination_sequences_simple", coordination_sequences::simple, (
      arg_("asu_mappings"),
      arg_("pair_asu_table_table"),
      arg_("n_shells")));
  }

}} // namespace coordination_sequences::<anonymous>

namespace boost_python {

  void
  wrap_coordination_sequences() { coordination_sequences::wrap_all(); }

}}} // namespace cctbx::crystal::boost_python
