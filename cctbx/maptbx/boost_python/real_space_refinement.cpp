// Done by Erik McKee
#include <cctbx/boost_python/flex_fwd.h>

#include <cctbx/maptbx/real_space_refinement.h>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/overloads.hpp>

namespace cctbx { namespace maptbx { namespace boost_python {

namespace {

  struct real_space_refinement_wrappers
  {
    BOOST_PYTHON_FUNCTION_OVERLOADS(
      real_space_refinement_gradients_overloads,
      real_space_refinement::gradients, 3, 5)

    static void
    wrap()
    {
      using namespace boost::python;
      def("real_space_refinement_residual",
        (double(*)(
          af::const_ref<double, af::flex_grid<> > const&,
          scitbx::mat3<double> const&,
          af::const_ref<scitbx::vec3<double> > const&,
          af::const_ref<double> const&))
        real_space_refinement::residual, (
        arg_("map"), arg_("gridding_matrix"), arg_("sites_cart"), arg_("weights")));
      def("real_space_refinement_gradients",
        (af::shared<scitbx::vec3<double> >(*)(
          af::const_ref<double, af::flex_grid<> > const&,
          scitbx::mat3<double> const&,
          af::const_ref<scitbx::vec3<double> > const&,
          double,std::size_t))
            real_space_refinement::gradients,
        real_space_refinement_gradients_overloads((
          arg_("map"), arg_("gridding_matrix"), arg_("sites_cart"),
          arg_("delta_h")=1.0,
          arg_("max_iter")=0)));
    }
  };

} // namespace <anoymous>

  void wrap_real_space_refinement()
  {
    real_space_refinement_wrappers::wrap();
  }

}}} // namespace cctbx::maptbx::boost_python
