#ifndef BUF_BASED_SVC_H
#define BUF_BASED_SVC_H
#include <cstddef>

namespace iotbx {
  namespace detectors {

void buffer_uncompress(char*, std::size_t, int*);

  }//namespace detectors
}//namespace iotbx

#endif
