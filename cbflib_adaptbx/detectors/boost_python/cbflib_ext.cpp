#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <cbflib_adaptbx/basic.h>
#include <cbflib_adaptbx/mar_adaptor.h>
#include <cbflib_adaptbx/cbf_adaptor.h>
#include <cbflib_adaptbx/sls_pilatus_adaptor.h>
#include <cbflib_adaptbx/general_cbf_write.h>

struct dummy {}; // work around gcc-3.3-darwin bug

#include <boost/python.hpp>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;
using namespace iotbx::detectors;

namespace iotbx{
namespace detectors{

bool assert_equal(scitbx::af::flex_int read1, scitbx::af::flex_int read2){
  SCITBX_ASSERT(read1.size()==read2.size());
  int N = read1.size();
  int* begin1 = read1.begin();
  int* begin2 = read2.begin();
  for (int n = 0; n<N; ++n){ SCITBX_ASSERT(*begin1++ == *begin2++); }
  return true;
}
}}

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
   class_<MiniCBFAdaptor, bases<CBFAdaptor> >("MiniCBFAdaptor",init<std::string>())
     .def("read_data",&MiniCBFAdaptor::read_data)
     .def("buffer_based_uncompress",&MiniCBFAdaptor::buffer_based_uncompress)
     .def("optimized_read_data",(scitbx::af::flex_int (MiniCBFAdaptor::*)(const int&, const int&))&MiniCBFAdaptor::optimized_read_data)
     .def("optimized_read_data",(scitbx::af::flex_int (MiniCBFAdaptor::*)())&MiniCBFAdaptor::optimized_read_data)
   ;
   class_<CBFWriteAdaptor, bases<CBFAdaptor> >("CBFWriteAdaptor",init<std::string>())
     .def("write_data",&CBFWriteAdaptor::write_data)
   ;
   def("assert_equal",&assert_equal);

}
