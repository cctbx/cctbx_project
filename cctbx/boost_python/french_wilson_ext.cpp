#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <boost/python/tuple.hpp>
#include <cctbx/french_wilson.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace boost_python {

  // (valid, I, SIGI, h, prior_dominated)
  boost::python::tuple
  invert_french_wilson_wrapper(
    double F, double SIGF, double mean_intensity, bool centric, double h_min)
  {
    double i_obs(0), sig_i_obs(0), h(0);
    bool prior_dominated(false);
    bool valid = invert_french_wilson(F, SIGF, mean_intensity, centric,
      i_obs, sig_i_obs, h, prior_dominated, h_min);
    return boost::python::make_tuple(
      valid, i_obs, sig_i_obs, h, prior_dominated);
  }

  void init_module()
  {
    using namespace boost::python;

    def("invert_french_wilson", invert_french_wilson_wrapper, (
      arg("F"),
      arg("SIGF"),
      arg("mean_intensity"),
      arg("centric"),
      arg("h_min")=-6.0));

    def("expectEFW",
      (double(*)
        (double,
         double,
         bool)) expectEFW, (
      arg("eosq"),
      arg("sigesq"),
      arg("centric")));

    def("expectEsqFW",
      (double(*)
        (double,
         double,
         bool)) expectEsqFW, (
      arg("eosq"),
      arg("sigesq"),
      arg("centric")));

    def("is_FrenchWilson",
      (bool(*)
        (af::shared<double>,
         af::shared<double>,
         af::shared<bool>,
         double)) is_FrenchWilson, (
      arg("F"),
      arg("SIGF"),
      arg("is_centric"),
      arg("eps")));

  }

}}

BOOST_PYTHON_MODULE(cctbx_french_wilson_ext)
{
  cctbx::boost_python::init_module();
}
