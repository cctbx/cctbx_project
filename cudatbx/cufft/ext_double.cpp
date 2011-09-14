
#ifndef CUFFT_DOUBLE_PRECSION
#define CUFFT_DOUBLE_PRECISION 1
#endif
#include <gputbx/cufft/cufft.hpp>

#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <scitbx/array_family/boost_python/utils.h>
#include <scitbx/array_family/versa.h>
#include <scitbx/array_family/accessors/c_grid.h>

namespace gputbx { namespace cufft {

  void wrap_cufft_double_precision ()
  {
    using namespace boost::python;
    def("real_to_complex_3d_in_place_dp", real_to_complex_3d_in_place, (
      arg("data")));
    def("complex_to_complex_3d_in_place_dp", complex_to_complex_3d_in_place, (
      arg("data"),
      arg("direction")));
    def("complex_to_real_3d_in_place_dp", complex_to_real_3d_in_place, (
      arg("data"),
      arg("n")));
  }

}} // namespace gputbx::cufft
