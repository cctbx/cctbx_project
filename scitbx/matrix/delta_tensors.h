#ifndef SCITBX_DELTA_TENSORS_H
#define SCITBX_DELTA_TENSORS_H

namespace scitbx { namespace matrix {

  template<typename ResultType, typename IntType>
  ResultType delta_x_delta(IntType i, IntType j, IntType k, IntType l) {
    if (i == j && k == l) return 1;
    return 0;
  }

}}

#endif
