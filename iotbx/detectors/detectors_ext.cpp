#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <scitbx/array_family/flex_types.h>

namespace af = scitbx::af;

namespace {

af::flex_int ReadADSC(const std::string& filename,
                      const long& ptr, const long& size1,
                      const long& size2,const int& big_endian ) {
  std::ifstream cin(filename.c_str());
  long fileLength = ptr + 2 * size1 * size2;
  std::vector<char> chardata(fileLength);
  cin.read(&*chardata.begin(),fileLength);
  cin.close();

  unsigned char* uchardata = (unsigned char*) &*chardata.begin();

  af::flex_int z(af::flex_grid<>(size1,size2));

  int* begin = z.begin();
  std::size_t sz = z.size();

  //af::ref<int> r(z.begin(),z.size());
  //the redesign; use r[i] = and r.size()

  if (big_endian) {
    for (std::size_t i = 0; i < sz; i++) {
      begin[i] = 256 * uchardata[ptr+2*i] + uchardata[ptr + 2*i +1];
    }
  } else {
    for (std::size_t i = 0; i < sz; i++) {
      begin[i] = 256 * uchardata[ptr+2*i+1] + uchardata[ptr + 2*i];
    }
  }

  return z;
}

af::flex_int ReadMAR(const std::string& filename,
                      const long& ptr, const long& size1,
                      const long& size2,const int& big_endian ) {
  return ReadADSC(filename,ptr,size1,size2,big_endian);
}

af::flex_int Bin2_by_2(const af::flex_int& olddata) {
  int oldsize = olddata.size();
  int olddim = std::sqrt((double)oldsize); // implicit real-valued intermediate
  int newdim = olddim/2;
  // always assume a square image!!!
  af::flex_int newdata(af::flex_grid<>(newdim,newdim));
  int *newptr = newdata.begin();

  const int *old = olddata.begin();
  for (std::size_t i = 0; i<newdim; ++i) { //think row-down
    for (std::size_t j = 0; j<newdim; ++j) { // think column-across
      *newptr++ =  old[2*i*olddim+2*j] + old[2*i*olddim+2*j+1] +
                         old[(2*i+1)*olddim+2*j] + old[(2*i+1)*olddim+2*j+1] ;
    }
  }
  return newdata;
}

struct dummy {}; // work around gcc-3.3-darwin bug

} // namespace <anonymous>

#include <boost/python.hpp>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;

BOOST_PYTHON_MODULE(iotbx_detectors_ext)
{
#if defined(__APPLE__) && defined(__MACH__) \
 && defined(__GNUC__) && __GNUC__ == 3 && __GNUC_MINOR__ == 3
   class_<dummy>("_dummy", no_init);
#endif
   def("ReadADSC", ReadADSC);
   def("ReadMAR", ReadMAR);
   def("Bin2_by_2", Bin2_by_2);
}
