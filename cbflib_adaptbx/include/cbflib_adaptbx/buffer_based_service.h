#ifndef BUF_BASED_SVC_H
#define BUF_BASED_SVC_H
#include <cstddef>
#include <vector>

namespace iotbx {
  namespace detectors {

void buffer_uncompress(const char*, std::size_t, int*);
std::vector<char> buffer_compress(const int*, const std::size_t&);

  }//namespace detectors
}//namespace iotbx

#endif
