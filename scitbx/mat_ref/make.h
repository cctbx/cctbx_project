#ifndef SCITBX_MAT_REF_MAKE_H
#define SCITBX_MAT_REF_MAKE_H

#include <scitbx/array_family/versa.h>
#include <scitbx/mat_ref.h>

namespace scitbx {

  template <typename NumType, class AccessorType>
  mat_const_ref<NumType, AccessorType>
  mat_const_ref_to(af::versa<NumType, AccessorType> const &a) {
    return mat_const_ref<NumType, AccessorType>(a.const_ref());
  }

  template <typename NumType, class AccessorType>
  mat_ref<NumType, AccessorType>
  mat_ref_to(af::versa<NumType, AccessorType> &a) {
    return mat_ref<NumType, AccessorType>(a.ref());
  }

} // scitbx

#endif // GUARD
