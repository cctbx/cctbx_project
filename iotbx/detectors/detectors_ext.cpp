#include <scitbx/array_family/boost_python/flex_fwd.h>
#include <string>
#include <vector>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <exception>
#include <scitbx/array_family/flex_types.h>
#include <scitbx/math/utils.h>

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

af::flex_int ReadRAXIS(const std::string& characters,
                       const int& width, const long& size1,
                       const long& size2,
                       const int& big_endian ) {
  af::flex_int z(af::flex_grid<>(size1,size2));

  int* begin = z.begin();
  std::size_t sz = z.size();

  std::string::const_iterator ptr = characters.begin();
  char* raw = new char[2];

  if (big_endian) {
    for (std::size_t i = 0; i < sz; i++) {
      raw[0] = *(ptr++); raw[1] = *(ptr++);
      unsigned short int* usi_raw = reinterpret_cast<unsigned short int*>(raw);
      if (*usi_raw <= 32767) {
        begin[i] = *usi_raw;
      } else {
        begin[i] = ((signed short int)(*usi_raw) + 32768) * 32;
      }
    }

  } else {
    for (std::size_t i = 0; i < sz; i++) {
      raw[1] = *(ptr++); raw[0] = *(ptr++);
      unsigned short int* usi_raw = reinterpret_cast<unsigned short int*>(raw);
      if (*usi_raw <= 32767) {
        begin[i] = *usi_raw;
      } else {
        begin[i] = ((signed short int)(*usi_raw) + 32768) * 32;
      }
    }
  }

  delete[] raw;
  return z;
}

std::string
unpad_raxis(const std::string& in, const int& recordlength, const int& pad){
  int number_records = in.size()/recordlength;
  int outlength = (recordlength-pad)*number_records;
  std::vector<char> out;
  out.reserve(outlength);
  std::string::const_iterator inptr=in.begin();
  for (int irec=0; irec<number_records; ++irec){
    for (int y=0; y<(recordlength-pad); ++y,++inptr) {
      out.push_back(*inptr);
    }
    inptr+=pad;
  }
  return std::string(out.begin(),out.end());
}

af::flex_int MakeSquareRAXIS(const int& np,
                             const int& extra,
                             const int& oldnslow,
                             const af::flex_int& rawlinearintdata) {
  /* In application, this function resizes the Raxis II image of 1900x1900
     rectangular pixels into a new image with 1962x1962 square pixels. The
     fast dimension is padded on each side with zeroes, while the slow 
     dimension is resampled so the data are stretched into the new size. 
  */
  af::flex_int z(af::flex_grid<>(np,np));

  int* begin = z.begin();
  const int* raw   = rawlinearintdata.begin();

  for (std::size_t slow = 0; slow < 2*extra; ++slow){
    for (std::size_t fast = 0; fast < np; ++fast){
      *(begin+slow*np+fast) = fast;}}
  for (std::size_t slow = 0; slow < oldnslow; ++slow){
    for (std::size_t fast = 0; fast < extra; ++fast){
      *(begin+(2*extra+slow)*np+fast) = 0;}
    for (std::size_t fast = 0; fast < oldnslow; ++fast){
      *(begin+(2*extra+slow)*np+extra+fast) = *(raw+(slow*oldnslow)+fast);}
    for (std::size_t fast = 0; fast < extra; ++fast){
      *(begin+(2*extra+slow)*np+extra+oldnslow+fast) = 0;}
  }

  //rearrange data in place
  for (std::size_t newslow = 0; newslow < np; ++newslow){
    double float_begin = (static_cast<double>(newslow)/np)*
                          static_cast<double>(oldnslow)+(2*extra);
    double float_end = (static_cast<double>(newslow+1)/np)*
                        static_cast<double>(oldnslow)+(2*extra);
    std::vector<int> signatures_addr;
    std::vector<double> signatures_frac;

    for (std::size_t xf = int(float_begin); xf < 1 + float_end; ++xf){
        double xf_min;
        if (float_begin < static_cast<double>(xf)){ xf_min=xf; }
        else { xf_min=float_begin; }

        double xf_max;
        if (float_end > static_cast<double>(xf+1)){xf_max=xf+1;}
        else { xf_max=float_end; }

        if (xf_max-xf_min>1.e-7) {
          signatures_addr.push_back(xf);
          signatures_frac.push_back(xf_max-xf_min);}
    }

    af::flex_int row_temp(af::flex_grid<>(np,1));
    int* rt = row_temp.begin();
    for (std::size_t item = 0; item < signatures_addr.size(); ++item){
      for (std::size_t fast = 0; fast < np; ++fast){
          *(rt+fast)+=
            *(begin + signatures_addr[item]*np + fast) * signatures_frac[item];
      }
      for (std::size_t fast = 0; fast < np; ++fast){
        *(begin + newslow*np + fast)=*(rt+fast);
      }
    }
  }

  return z;
}

af::flex_int Bin2_by_2(const af::flex_int& olddata) {
  int oldsize = olddata.size();
  int olddim = std::sqrt((double)oldsize);
  int newdim = olddim/2;
  if (olddim%2!=0) {throw;} // image dimension must be even so it can be divided by 2
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
   def("ReadRAXIS", ReadRAXIS);
   def("unpad_raxis", unpad_raxis);
   def("MakeSquareRAXIS", MakeSquareRAXIS);
   def("Bin2_by_2", Bin2_by_2);
}
