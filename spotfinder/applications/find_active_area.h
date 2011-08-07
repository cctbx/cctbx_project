#ifndef FIND_ACTIVE_AREA_H
#define FIND_ACTIVE_AREA_H

#include <cstring>
#include <cmath>

#include <scitbx/error.h>//SCITBX_CHECK_POINT
#include <scitbx/array_family/shared.h>
#include <scitbx/array_family/flex_types.h>

namespace af = scitbx::af;
namespace spotfinder {
namespace distltbx {

scitbx::af::shared<int>
find_active_area (scitbx::af::flex_int const& rawdata) {
  std::size_t size1(rawdata.accessor().focus()[0]); //slow
  std::size_t size2(rawdata.accessor().focus()[1]); //fast
  af::flex_int z(scitbx::af::flex_grid<>(size1+2,size2+2));
  const int* raw = rawdata.begin();
  int* pad = z.begin();

  //copy the raw data into a temporary array padded on all edges by zero
  for (size_t islow=0; islow<size1; ++islow){
     const void* src = (const void*)( raw + islow*size2 );
     void* dst = (void*)( pad + (islow +1)*(size2+2) + 1 );
     std::memcpy( dst, src, size2*sizeof(int) );
  }
  scitbx::af::shared<int> result;

  size_t pad_data_begin = size2+3;
  size_t pad_data_end = (size2+2)*(size1+1);
  // test upper left corner of active area

  for (const int* dst = pad + pad_data_begin; dst < pad + pad_data_end; ++dst){
    if (*dst!=0 && *(dst-1)==0 && *(dst+size2+1)==0 &&
        *(dst-size2-2)==0 && *(dst-size2-3)==0){
      //decode unpadded location
      int slow = (((dst - pad)/(size2+2))-1 );
      int fast = (((dst - pad)%(size2+2))-1 );
      result.push_back( slow );
      result.push_back( fast );
    }
  }

  // test lower right corner of active area

  for (const int* dst = pad + pad_data_begin; dst < pad + pad_data_end; ++dst){
    if (*dst!=0 && *(dst+1)==0 && *(dst-size2-1)==0 &&
        *(dst+size2+1)==0 && *(dst+size2+2)==0 && *(dst+size2+3)==0){
      //decode unpadded location
      int slow = -(((dst - pad)/(size2+2))-1 );
      int fast = -(((dst - pad)%(size2+2))-1 );
      result.push_back( slow );//negative number signifies lower right
      result.push_back( fast );
    }
  }

  return result;
};

} //namespace

} //namespace

#endif //FIND_ACTIVE_AREA_H
