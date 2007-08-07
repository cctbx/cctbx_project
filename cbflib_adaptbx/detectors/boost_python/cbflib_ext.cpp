#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cbflib_adaptbx/basic.h>
#include <cbflib_adaptbx/mar_adaptor.h>
#include <cbflib_adaptbx/cbf_adaptor.h>

struct dummy {}; // work around gcc-3.3-darwin bug

#include <boost/python.hpp>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;
using namespace iotbx::detectors;

BOOST_PYTHON_MODULE(cbflib_ext)
{
#if defined(__APPLE__) && defined(__MACH__) \
 && defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 3
   class_<dummy>("_dummy", no_init);
#endif
   class_<Mar345Adaptor >("Mar345Adaptor",init<std::string>())
     .def("read_header",&Mar345Adaptor::read_header)
     .def("read_data",&Mar345Adaptor::read_data)
     .def("pixel_size",&Mar345Adaptor::pixel_size)
     .def("wavelength",&Mar345Adaptor::wavelength)
     .def("distance",&Mar345Adaptor::distance)
     .def("gain",&Mar345Adaptor::gain)
     .def("overload",&Mar345Adaptor::overload)
     .def("size1",&Mar345Adaptor::size1)
     .def("size2",&Mar345Adaptor::size2)
     .def("osc_range",&Mar345Adaptor::osc_range)
     .def("osc_start",&Mar345Adaptor::osc_start)
     .def("twotheta",&Mar345Adaptor::twotheta)
     .def("rawdata",&Mar345Adaptor::rawdata)
     .def("test",&Mar345Adaptor::test)
   ;
   class_<CBFAdaptor >("CBFAdaptor",init<std::string>())
     .def("read_header",&CBFAdaptor::read_header)
     .def("read_data",&CBFAdaptor::read_data)
     .def("rawdata",&CBFAdaptor::read_data)
     .def("pixel_size",&CBFAdaptor::pixel_size)
     .def("wavelength",&CBFAdaptor::wavelength)
     .def("distance",&CBFAdaptor::distance)
     //.def("gain",&CBFAdaptor::gain)
     .def("overload",&CBFAdaptor::overload)
     .def("size1",&CBFAdaptor::size1)
     .def("size2",&CBFAdaptor::size2)
     .def_readonly("beam_index_slow",&CBFAdaptor::beam_index1)
     .def_readonly("beam_index_fast",&CBFAdaptor::beam_index2)
     .def("osc_range",&CBFAdaptor::osc_range)
     .def("osc_start",&CBFAdaptor::osc_start)
     .def("twotheta",&CBFAdaptor::twotheta)
     .def("raster_description",&CBFAdaptor::raster_description)
     //.def("test",&CBFAdaptor::test)
   ;
}
