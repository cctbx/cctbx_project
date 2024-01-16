#include "kokkostbx/kokkos_utils.h"

namespace kokkostbx {
void transfer_double2kokkos(vector_cudareal_t& dst, const double* src, const size_t length) {
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

void transfer_vector2kokkos(view_1d_t<bool>& dst, const std::vector<bool>& src) {
    if (true) {
        // printf("== Transfer %s from %p\n", dst.label().c_str(), (void*) dst.data());
        // printf(" - size src|dst: %d|%d\n", src.size(), dst.span() );
    }
    if (dst.span() < src.size()) {
        resize(dst, src.size());
        // printf(" - size changed, new size: %d\n", dst.span() );
    }
    auto host_view = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < src.size(); ++i) {
        host_view(i) = src[i];
    }
    Kokkos::deep_copy(dst, host_view);
}

}  // namespace kokkostbx
