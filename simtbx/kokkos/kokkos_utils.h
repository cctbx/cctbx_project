#ifndef SIMTBX_KOKKOS_UTILS_H
#define SIMTBX_KOKKOS_UTILS_H

#include "scitbx/array_family/shared.h"
#include "scitbx/array_family/flex_types.h"
#include "kokkos_types.h"

namespace {
  template <typename T, typename U>
  void
  transfer_2kokkos(view_1d_t<T> &dst, const U &src) {
    auto host_view = Kokkos::create_mirror_view(dst);
    auto src_ptr = src.begin();
    for (int i=0; i<src.size(); ++i) {
      host_view( i ) = src_ptr[ i ];
    }
    Kokkos::deep_copy(dst, host_view);
  }

  template <typename T, typename U>
  void
  transfer_kokkos2(T &dst, const view_1d_t<U> &src) {
    auto host_view = Kokkos::create_mirror_view(src);
    Kokkos::deep_copy(host_view, src);
    auto dst_ptr = dst.begin();
    for (int i=0; i<host_view.span(); ++i) {
      dst_ptr[ i ] = host_view( i );
    }
  }
}

namespace simtbx { namespace Kokkos { 

namespace af = scitbx::af;

  template <typename T>
  void
  transfer_flex2kokkos(view_1d_t<T> &dst, const af::shared<T>  &src) {
    if (true) {
      printf("== Transfer %s from %p\n", dst.label().c_str(), (void*) dst.data());
      printf(" - size src|dst: %d|%d\n", src.size(), dst.span() );
    }
    if (dst.span() < src.size()) {
      resize(dst, src.size());
      printf(" - size changed, new size: %d\n", dst.span() );
    }
    auto host_view = create_mirror_view(dst);

    for (int i=0; i<src.size(); ++i) {
      host_view( i ) = src[ i ];
    }
    deep_copy(dst, host_view);
    printf("copied from %p to %p.\n", (void*) host_view.data(), (void*) dst.data());

    print_view(dst);
  }

  template <typename T>
  void
  transfer_kokkos2flex(af::shared<T> &dst, const view_1d_t<T> &src) {
    if (true) {
      printf("== Transfer %s from %p\n", src.label().c_str(), (void*) src.data());
    }

    auto host_view = create_mirror_view(src);
    deep_copy(host_view, src);

    auto dst_ptr = dst.begin();
    for (int i=0; i<src.span(); ++i) {
      dst_ptr[ i ] = host_view( i );
    }
  }

} // Kokkos
} // simtbx

#endif // SIMTBX_KOKKOS_UTILS_H