#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_CUDA
#define MemSpace Kokkos::CudaSpace
#endif
#ifdef KOKKOS_ENABLE_HIP
#define MemSpace Kokkos::Experimental::HIPSpace
#endif
#ifdef KOKKOS_ENABLE_OPENMPTARGET
#define MemSpace Kokkos::OpenMPTargetSpace
#endif

#ifndef MemSpace
#define MemSpace Kokkos::HostSpace
#endif

using ExecSpace = MemSpace::execution_space;
using range_policy = Kokkos::RangePolicy<ExecSpace>;

// Allocate y, x vectors and Matrix A on device.
typedef Kokkos::View<double*, Kokkos::LayoutLeft, MemSpace>   ViewVectorType;
typedef Kokkos::View<double**, Kokkos::LayoutLeft, MemSpace>  ViewMatrixType;

template<typename T>
using view_1d_t = Kokkos::View<T*, MemSpace>;

using vector_bool_t = view_1d_t<bool>;
using vector_double_t = view_1d_t<double>;
using vector_float_t = view_1d_t<float>;
using vector_cudareal_t = view_1d_t<CUDAREAL>;
using vector_int_t = view_1d_t<int>;
using vector_ushort_t = view_1d_t<unsigned int short>;

//ViewVectorType y( "y", N );
//ViewVectorType x( "x", M );
