#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/miller/sym_equiv.h>
#include <cctbx/miller/math.h>
#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/class.hpp>
#include <scitbx/boost_python/container_conversions.h>

namespace cctbx { namespace miller { namespace boost_python {

  void wrap_asu();
  void wrap_bins();
  void wrap_expand_to_p1();
  void wrap_index_generator();
  void wrap_index_span();
  void wrap_match_bijvoet_mates();
  void wrap_match_indices();
  void wrap_merge_equivalents();
  void wrap_phase_transfer();
  void wrap_sym_equiv();

namespace {

  void register_tuple_mappings()
  {
    using namespace scitbx::boost_python::container_conversions;

    tuple_mapping_variable_capacity<af::shared<sym_equiv_index> >();
  }

  void init_module()
  {
    using namespace boost::python;

    register_tuple_mappings();

    wrap_sym_equiv(); // must be wrapped first to enable use of bases<>
    wrap_asu();
    wrap_bins();
    wrap_expand_to_p1();
    wrap_index_generator();
    wrap_index_span();
    wrap_match_bijvoet_mates();
    wrap_match_indices();
    wrap_merge_equivalents();
    wrap_phase_transfer();

    def("statistical_mean",
      (double(*)(sgtbx::space_group const&,
                 bool,
                 af::const_ref<index<> > const&,
                 af::const_ref<double> const&)) statistical_mean);
  }

} // namespace <anonymous>
}}} // namespace cctbx::miller::boost_python

BOOST_PYTHON_MODULE(cctbx_miller_ext)
{
  cctbx::miller::boost_python::init_module();
}
