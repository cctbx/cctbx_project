#ifndef KOKKOS_UTILS_H
#define KOKKOS_UTILS_H

#include <vector>
#include <algorithm>

#include "kokkostbx/kokkos_types.h"
#include "scitbx/array_family/flex_types.h"
#include "scitbx/array_family/shared.h"

//*******************************************************************
// Transfer data
//*******************************************************************
namespace {
template <typename T, typename U>
void transfer_X2kokkos(view_1d_t<T>& dst, const U& src) {
    auto host_view = Kokkos::create_mirror_view(dst);
    auto src_ptr = src.begin();
    for (int i = 0; i < src.size(); ++i) {
        host_view(i) = src_ptr[i];
    }
    Kokkos::deep_copy(dst, host_view);
}

template <typename T, typename U>
void transfer_X2kokkos(view_1d_t<T>& dst, const U* src, const size_t length) {
    auto host_view = Kokkos::create_mirror_view(dst);
    for (int i = 0; i < length; ++i) {
        host_view(i) = src[i];
    }
    Kokkos::deep_copy(dst, host_view);
}

template <typename T, typename U>
void transfer_kokkos2X(T& dst, const view_1d_t<U>& src) {
    auto host_view = Kokkos::create_mirror_view(src);
    Kokkos::deep_copy(host_view, src);
    auto dst_ptr = dst.begin();
    for (int i = 0; i < host_view.span(); ++i) {
        dst_ptr[i] = host_view(i);
    }
}
}  // namespace

//*******************************************************************
// math functions
//*******************************************************************

namespace kokkostbx {

namespace af = scitbx::af;

void transfer_double2kokkos(vector_cudareal_t& dst, const double* src, const size_t length);

template <typename T>
void transfer_shared2kokkos(view_1d_t<T>& dst, const af::shared<T>& src) {
    if (true) {
        // printf("== Transfer %s from %p\n", dst.label().c_str(), (void*) dst.data());
        // printf(" - size src|dst: %d|%d\n", src.size(), dst.span() );
    }
    if (dst.span() < src.size()) {
        resize(dst, src.size());
        // printf(" - size changed, new size: %d\n", dst.span() );
    }

    transfer_X2kokkos(dst, src);
    if (src.size() < 20) {
        // print_view(dst);
    }
}

template <typename T>
void transfer_kokkos2shared(af::shared<T>& dst, const view_1d_t<T>& src) {
    if (true) {
        // printf("== Transfer %s from %p\n", src.label().c_str(), (void*)
        // src.data());
    }
    transfer_kokkos2X(dst, src);
}

void transfer_vector2kokkos(view_1d_t<bool>& dst, const std::vector<bool>& src);

template <typename T>
void transfer_vector2kokkos(view_1d_t<T>& dst, const std::vector<T>& src) {
    if (true) {
        // printf("== Transfer %s from %p\n", dst.label().c_str(), (void*) dst.data());
        // printf(" - size src|dst: %d|%d\n", src.size(), dst.span() );
    }
    if (dst.span() < src.size()) {
        resize(dst, src.size());
        // printf(" - size changed, new size: %d\n", dst.span() );
    }
    auto host_view = Kokkos::View<const T*, Kokkos::HostSpace>(src.data(), src.size());
    auto dst_subview = Kokkos::subview(dst, std::pair<int, int>(0, src.size()));
    Kokkos::deep_copy(dst_subview, host_view);
}

template <typename T>
void transfer_kokkos2vector(std::vector<T>& dst, const view_1d_t<T>& src) {
    auto host_view = Kokkos::create_mirror_view(src);
    // printf(" - size src|dst: %d|%d\n", src.span(), dst.size() );
    const int length = std::min(dst.size(), src.span());
    Kokkos::deep_copy(host_view, src);
    for (int i = 0; i < length; ++i) {
        dst[i] = host_view(i);
    }
}

template <typename T>
void transfer_kokkos2flex(af::flex_double& dst, const view_1d_t<T>& src) {
    if (true) {
        // printf("== Transfer %s from %p\n", src.label().c_str(), (void*)
        // src.data());
    }
    transfer_kokkos2X(dst, src);
}

template <typename T>
void transfer_kokkos2flex(af::flex_int& dst, const view_1d_t<T>& src) {
    if (true) {
        // printf("== Transfer %s from %p\n", src.label().c_str(), (void*)
        // src.data());
    }
    transfer_kokkos2X(dst, src);
}

}  // namespace kokkostbx

#endif  // KOKKOS_UTILS_H
