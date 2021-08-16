#include "simtbx/kokkos/kokkos_utils.h"

namespace simtbx { namespace Kokkos {
  void
  transfer_double2kokkos(vector_cudareal_t &dst, const double *src, const size_t length) {
    if (true) {
      // printf("== Transfer %s from %p\n", dst.label().c_str(), (void*) dst.data());
      // printf(" - size src|dst: %d|%d\n", length, dst.span() );
    }
    if (dst.span() < length) {
      resize(dst, length);
      // printf(" - size changed, new size: %d\n", dst.span() );
    }

    transfer_X2kokkos(dst, src, length);

    if (length < 20) {
      // print_view(dst);
    }
  }

} // Kokkos
} // simtbx
