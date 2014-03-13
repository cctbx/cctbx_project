#include <boost/python/module.hpp>
#include <boost/python/def.hpp>
#include <cctbx/french_wilson.h>
#include <scitbx/array_family/boost_python/shared_wrapper.h>
#include <cctbx/import_scitbx_af.h>

namespace cctbx { namespace boost_python {

  void init_module()
  {
    using namespace boost::python;

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

BOOST_PYTHON_MODULE(cctbx_ext)
{
  cctbx::boost_python::init_module();
}
