#ifndef CCTBX_XRAY_THICKNESS_H
#define CCTBX_XRAY_THICKNESS_H

#include <scitbx/array_family/shared.h>

namespace cctbx { namespace xray {

  template <typename FloatType>
  struct thickness {
    FloatType value;
    bool grad;
    int grad_index;
    /* a generic attribute that allows to pass information about this particular
    fraction. In the case of ED will keep index to the associated zone/frame
    index
    */
    int tag;

    thickness(FloatType value, bool grad=false)
      : value(value), grad(grad), grad_index(-1), tag(-1)
    {}
    thickness(FloatType value, int tag, bool grad=false)
      : value(value), grad(grad), grad_index(-1), tag(tag)
    {}
    thickness deep_copy() { return thickness(*this); }
  };

}} // namespace cctbx::xray

#endif // GUARD
