#ifndef CCTBX_XRAY_PACKING_ORDER_H
#define CCTBX_XRAY_PACKING_ORDER_H

namespace cctbx { namespace xray {

  template <int Version>
  struct packing_order_convention {};

  template <>
  struct packing_order_convention<1>
  {
    typedef bool check_version_at_compile_time;
  };

}} // namespace cctbx::xray

#endif // CCTBX_XRAY_PACKING_ORDER_H
