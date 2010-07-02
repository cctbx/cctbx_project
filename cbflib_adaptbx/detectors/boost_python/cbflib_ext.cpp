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

boost::python::str compressed_string(MiniCBFAdaptor& ada){
  //later, reorganize MiniCBFAdaptor so some of this code is in the class proper.
  ada.common_file_access();

  //had to make ada.cbf_h public for this
  wrapper_of_byte_decompression wrap_dee(&ada.cbf_h,ada.dim1*ada.dim2);
  wrap_dee.set_file_position();

  scitbx::af::shared<char> compressed_buffer(wrap_dee.size_text);
  char* buffer_begin = compressed_buffer.begin();
  std::size_t sz_buffer = compressed_buffer.size();

  wrap_dee.copy_raw_compressed_string_to_buffer(buffer_begin, sz_buffer);
  return boost::python::str(buffer_begin,sz_buffer);
}

scitbx::af::flex_int uncompress(const boost::python::str& packed,const int& slow, const int& fast){

    std::string strpacked = extract<std::string>(packed);
    std::size_t sz_buffer = strpacked.size();

    //C++ weirdness
    scitbx::af::flex_int z((scitbx::af::flex_grid<>(slow,fast)),scitbx::af::init_functor_null<int>());
    int* begin = z.begin();
    std::size_t sz = z.size();

    iotbx::detectors::buffer_uncompress(strpacked.c_str(), sz_buffer, begin);

    return z;
}

boost::python::str compress(const scitbx::af::flex_int z){

    const int* begin = z.begin();
    std::size_t sz = z.size();

    std::vector<char> packed = iotbx::detectors::buffer_compress(begin, sz);

    return boost::python::str(&*packed.begin(),packed.size());
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
     .def("optimized_read_data",
     (scitbx::af::flex_int (MiniCBFAdaptor::*)(const int&, const int&))&MiniCBFAdaptor::optimized_read_data)
     .def("optimized_read_data",
     (scitbx::af::flex_int (MiniCBFAdaptor::*)())&MiniCBFAdaptor::optimized_read_data)
     .def("compressed_string",&compressed_string)
     .def("dim_slow",&MiniCBFAdaptor::dim_slow)
     .def("dim_fast",&MiniCBFAdaptor::dim_fast)
   ;
   class_<CBFWriteAdaptor, bases<CBFAdaptor> >("CBFWriteAdaptor",init<std::string>())
     .def("write_data",&CBFWriteAdaptor::write_data)
   ;
   def("assert_equal",&assert_equal);
   def("uncompress",&uncompress,(
     arg_("packed"),
     arg_("slow"),
     arg_("fast")
     )
   );
   def("compress",&compress);

}
