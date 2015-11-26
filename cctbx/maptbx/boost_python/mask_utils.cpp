
#include <scitbx/array_family/accessors/flex_grid.h>
#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/mask_utils.h>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/copy_const_reference.hpp>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>

namespace cctbx { namespace maptbx {
  namespace {

  void wrap_sample_mask_regions_impl()
  {
    using namespace boost::python;
    // typedef return_value_policy<copy_const_reference> ccr;
    def("sample_mask_regions",
      (af::shared<scitbx::vec3<double> > (*)
        (af::const_ref<int, af::flex_grid<> > const&,
          int,
          int,
          int,
          cctbx::uctbx::unit_cell const&))
      sample_mask_regions, (
        arg("mask"),
        arg("n_zone"),
        arg("volume"),
        arg("sampling_rate"),
        arg("unit_cell")));
  }
} // namespace <anonymous>

namespace boost_python {

  void
  wrap_sample_mask_regions() { wrap_sample_mask_regions_impl(); }

}}} // namespace iotbx::pdb::boost_python
