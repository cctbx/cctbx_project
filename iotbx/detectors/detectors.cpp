#include <string>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/array_family/boost_python/flex_fwd.h>

namespace af = scitbx::af;

af::flex_int ReadADSC(const std::string& filename,
                      const long& ptr, const long& size1,
                      const long& size2) {
  std::ifstream cin(filename.c_str());
  long fileLength = ptr + 2 * size1 * size2;
  char* chardata = new char[fileLength];
  cin.read(chardata,fileLength);
  cin.close();

  unsigned char* uchardata = (unsigned char*) chardata;

  af::flex_int z(af::flex_grid<>(size1,size2));

  int* begin = z.begin();
  std::size_t sz = z.size();

  //af::ref<int> r(z.begin(),z.size());
  //the redesign; use r[i] = and r.size()

  for (std::size_t i = 0; i < sz; i++) {
    begin[i] = 256 * uchardata[ptr+2*i] + uchardata[ptr + 2*i +1];
  }

  delete[] chardata;
  return z;
}

#include <boost/python.hpp>
#include <scitbx/boost_python/utils.h>
using namespace boost::python;

BOOST_PYTHON_MODULE(detectors)
{
   //scitbx::boost_python::import_module(
   //   "scitbx_boost.array_family.flex_scitbx_ext");
   //import in the __init__ file instead
   def("ReadADSC", ReadADSC);
}
