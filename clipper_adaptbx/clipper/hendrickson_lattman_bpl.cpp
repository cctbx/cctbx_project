#include <cctbx/boost_python/flex_fwd.h>
#include <boost/python/class.hpp>
#include <boost/python/args.hpp>

#if defined(__linux__) && defined(__GNUC__) \
 && __GNUC__ == 2 && __GNUC_MINOR__ == 96
# undef isnan
#endif

#include <clipper/clipper.h>
#include <clipper/clipper-contrib.h>
#include <clipper/cctbx/clipper_cctbx.h>

namespace clipper { namespace {

  class Compute_phifom_from_abcd_interface
  {
    public:
      Compute_phifom_from_abcd_interface(
        cctbx::uctbx::unit_cell const& unit_cell,
        cctbx::sgtbx::space_group const& space_group,
        af::shared<cctbx::miller::index<> > const& miller_indices,
        af::const_ref<cctbx::hendrickson_lattman<> > const&
          phase_probabilities)
      :
        miller_indices_(miller_indices),
        hkls_(CCTBX::as_HKL_info(
          unit_cell, space_group, miller_indices.const_ref())),
        phiw_(hkls_)
      {
        CCTBX_ASSERT(hkls_.num_reflections() == miller_indices.size());
        HKL_data<data64::ABCD> abcd = CCTBX::as_HKL_data(
          hkls_, miller_indices_.const_ref(), phase_probabilities);
        phiw_.compute(abcd, data64::Compute_phifom_from_abcd());
      }

      af::shared<double>
      centroid_phases() const
      {
        return CCTBX::extract_centroid_phases(
          phiw_, miller_indices_.const_ref());
      }

      af::shared<double>
      figures_of_merit() const
      {
        return CCTBX::extract_figures_of_merit(
          phiw_, miller_indices_.const_ref());
      }

    protected:
      af::shared<cctbx::miller::index<> > miller_indices_;
      HKL_info hkls_;
      HKL_data<data64::Phi_fom> phiw_;
  };

  struct Compute_phifom_from_abcd_interface_wrappers
  {
    typedef Compute_phifom_from_abcd_interface w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("Compute_phifom_from_abcd_interface", no_init)
        .def(init<cctbx::uctbx::unit_cell const&,
                  cctbx::sgtbx::space_group const&,
                  af::shared<cctbx::miller::index<> > const&,
                  af::const_ref<cctbx::hendrickson_lattman<> > const&>((
          arg_("unit_cell"),
          arg_("space_group"),
          arg_("miller_indices"),
          arg_("phase_probabilities"))))
        .def("centroid_phases", &w_t::centroid_phases)
        .def("figures_of_merit", &w_t::figures_of_merit)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_hendrickson_lattman()
  {
    Compute_phifom_from_abcd_interface_wrappers::wrap();
  }

}} // namespace clipper::boost_python
