// The core of this code is based on Kevin Cowtan's
// clipper/examples/csigmaa.cpp

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

  class SFweight_spline_interface
  {
    public:
      SFweight_spline_interface(
        cctbx::uctbx::unit_cell const& unit_cell,
        cctbx::sgtbx::space_group const& space_group,
        af::shared<cctbx::miller::index<> > const& miller_indices,
        bool anomalous_flag,
        af::const_ref<double> const& f_obs_data,
        af::const_ref<double> const& f_obs_sigmas,
        af::const_ref<std::complex<double> > const& f_calc,
        af::const_ref<bool> const& test_set_flags,
        int n_refln=1000,
        int n_param=20)
      :
        miller_indices_(miller_indices),
        hkls_(CCTBX::as_HKL_info(
          unit_cell, space_group, miller_indices.const_ref())),
        fb_(hkls_),
        fd_(hkls_),
        phiw_(hkls_),
        sfweight_spline_(n_refln, n_param)
      {
        if (anomalous_flag) {
          throw cctbx::error(
            "SFweight_spline_interface cannot currently"
            " handle anomalous arrays.");
        }
        CCTBX_ASSERT(test_set_flags.size() == miller_indices.size());
        CCTBX_ASSERT(hkls_.num_reflections() == miller_indices.size());
        af::const_ref<cctbx::miller::index<> >
          mi_ref = miller_indices_.const_ref();
        HKL_data<data64::F_sigF> f_obs_ = CCTBX::as_HKL_data(
          hkls_, mi_ref, f_obs_data, f_obs_sigmas);
        HKL_data<data64::F_phi> f_calc_ = CCTBX::as_HKL_data(
          hkls_, mi_ref, f_calc);
        HKL_data<data64::Flag> flag(hkls_);
        for(std::size_t i=0;i<mi_ref.size();i++) {
          data64::Flag fl;
          for(std::size_t i=0;i<mi_ref.size();i++) {
            HKL h = CCTBX::Hkl(mi_ref[i]);
            fl.flag() = test_set_flags[i] ? SFweight_spline<double>::BOTH
                                          : SFweight_spline<double>::NONE;
            flag.set_data(h, fl);
          }
        }
        // do sigmaa calc
        CCTBX_ASSERT(sfweight_spline_(fb_, fd_, phiw_, f_obs_, f_calc_, flag));
      }

      std::size_t
      number_of_spline_parameters()
      {
        return sfweight_spline_.params_s().size();
      }

      af::shared<std::complex<double> >
      fb() const
      {
        return CCTBX::extract_complex(fb_, miller_indices_.const_ref());
      }

      af::shared<std::complex<double> >
      fd() const
      {
        return CCTBX::extract_complex(fd_, miller_indices_.const_ref());
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
      HKL_data<data64::F_phi> fb_;
      HKL_data<data64::F_phi> fd_;
      HKL_data<data64::Phi_fom> phiw_;
      SFweight_spline<double> sfweight_spline_;
  };

  struct SFweight_spline_interface_wrappers
  {
    typedef SFweight_spline_interface w_t;

    static void
    wrap()
    {
      using namespace boost::python;
      class_<w_t>("SFweight_spline_interface", no_init)
        .def(init<cctbx::uctbx::unit_cell const&,
                  cctbx::sgtbx::space_group const&,
                  af::shared<cctbx::miller::index<> > const&,
                  bool,
                  af::const_ref<double> const&,
                  af::const_ref<double> const&,
                  af::const_ref<std::complex<double> > const&,
                  af::const_ref<bool> const&,
                  optional<int, int> >((
          arg_("unit_cell"),
          arg_("space_group"),
          arg_("miller_indices"),
          arg_("anomalous_flag"),
          arg_("f_obs_data"),
          arg_("f_obs_sigmas"),
          arg_("f_calc"),
          arg_("test_set_flags"),
          arg_("n_refln")=1000,
          arg_("n_param")=20)))
        .def("number_of_spline_parameters", &w_t::number_of_spline_parameters)
        .def("fb", &w_t::fb)
        .def("fd", &w_t::fd)
        .def("centroid_phases", &w_t::centroid_phases)
        .def("figures_of_merit", &w_t::figures_of_merit)
      ;
    }
  };

} // namespace <anonymous>

namespace boost_python {

  void wrap_sigmaa()
  {
    SFweight_spline_interface_wrappers::wrap();
  }

}} // namespace clipper::boost_python
