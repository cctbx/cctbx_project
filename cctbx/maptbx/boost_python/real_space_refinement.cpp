// Done by Erik McKee
#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/real_space_refinement.h>
#include <boost/python/def.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct real_space_refinement_wrappers
  {
    static void
    wrap()
    {
      using namespace boost::python;
      def("real_space_refinement_residual", real_space_refinement::residual<double>, (
        arg_("map"), arg_("gridding_matrix"), arg_("sites_cart")));
      def("real_space_refinement_gradients", real_space_refinement::gradients<double>, (
        arg_("map"), arg_("gridding_matrix"), arg_("sites_cart")));
    }
  };

} // namespace <anoymous>

  void wrap_real_space_refinement()
  {
    real_space_refinement_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
