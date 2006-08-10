#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <exception>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/math/utils.h>
#include <iotbx/detectors/bruker.h>

namespace af = scitbx::af;

namespace {
struct dummy {}; // work around gcc-3.3-darwin bug
} // namespace <anonymous>

#include <boost/python.hpp>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;

BOOST_PYTHON_MODULE(iotbx_detectors_bruker_ext)
{
#if defined(__APPLE__) && defined(__MACH__) \
 && defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 3
   class_<dummy>("_dummy", no_init);
#endif
   class_<iotbx::detectors::bruker>("Bruker_base", init<std::string>() )
   .def("linearintdata",&iotbx::detectors::bruker::linearintdata)
   .def_readonly("ccd_image_saturation",&iotbx::detectors::bruker::saturate)
   .def_readonly("pixel_size",&iotbx::detectors::bruker::pixsizemm)
   .def_readonly("osc_start",&iotbx::detectors::bruker::osc_start)
   .def_readonly("distance_cm",&iotbx::detectors::bruker::distance)
   .def_readonly("distance_delta",&iotbx::detectors::bruker::delta)
   .def_readonly("wavelength",&iotbx::detectors::bruker::wavelen)
   .def_readonly("centerx_pix",&iotbx::detectors::bruker::centerx)
   .def_readonly("centery_pix",&iotbx::detectors::bruker::centery)
   .def_readonly("osc_range",&iotbx::detectors::bruker::oscrange)
   .def_readonly("twotheta",&iotbx::detectors::bruker::twoth)
   ;
}
