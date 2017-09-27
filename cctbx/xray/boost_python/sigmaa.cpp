#include <cctbx/boost_python/flex_fwd.h>

#include <boost/python/class.hpp>
#include <boost/python/def.hpp>
#include <boost/python/args.hpp>
#include <boost/python/return_value_policy.hpp>
#include <boost/python/return_by_value.hpp>
#include <cctbx/xray/sigmaa.h>

namespace cctbx { namespace xray { namespace boost_python {

  void wrap_sigmaa()
  {
    using namespace boost::python;

    def("compute",
      (af::shared<double>(*)(
        af::const_ref<double> const&,
        af::const_ref<std::complex<double> > const&)) sigmaa::compute, (
          arg("f_obs"),
          arg("f_calc")));

  }

}}} // namespace cctbx::xray::boost_python
