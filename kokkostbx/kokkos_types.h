#ifndef KOKKOS_TYPES_H
#define KOKKOS_TYPES_H
#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_CUDA
    #define MemSpace Kokkos::CudaSpace
#endif
#ifdef KOKKOS_ENABLE_HIP
    #define MemSpace Kokkos::HIPSpace
#endif
#ifdef KOKKOS_ENABLE_OPENMPTARGET
    #define MemSpace Kokkos::OpenMPTargetSpace
#endif

#ifndef MemSpace
    #define MemSpace Kokkos::HostSpace
#endif

using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy<ExecSpace>;

template <typename T> using view_1d_t = Kokkos::View<T*, MemSpace>;
template <typename T> using view_4d_t = Kokkos::View<T****, MemSpace>;
template <typename T> using view_5d_t = Kokkos::View<T*****, MemSpace>;
template <typename T> using view_6d_t = Kokkos::View<T******, MemSpace>;
template <typename T> using view_6d6_t = Kokkos::View<T******[6], MemSpace>;

using vector_bool_t = view_1d_t<bool>;
using vector_double_t = view_1d_t<double>;
using vector_float_t = view_1d_t<float>;
using vector_cudareal_t = view_1d_t<CUDAREAL>;
using vector_int_t = view_1d_t<int>;
using vector_size_t = view_1d_t<std::size_t>;
using vector_uint_t = view_1d_t<unsigned int>;
using vector_ushort_t = view_1d_t<unsigned int short>;

template <typename T>
void print_view(const view_1d_t<T>& arg_view, size_t arg_first, size_t arg_last) {
    std::string label = arg_view.label();
    printf("print_view '%s'\n", label.c_str());

    Kokkos::parallel_for(
        "print_view", range_policy(arg_first, arg_last),
        KOKKOS_LAMBDA(const int i) { printf(" P[ %d ] = '%f'\n", i, (double)arg_view(i)); });
    Kokkos::fence();
}

template <typename T>
void print_view(const view_1d_t<T>& arg_view) {
    print_view(arg_view, 0, arg_view.span());
}

#endif  // KOKKOS_TYPES_H
