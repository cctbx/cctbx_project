#ifndef KOKKOSTBX_VECTOR2_H
#define KOKKOSTBX_VECTOR2_H

#include "kokkos_vector.h"

// #ifndef KOKKOS_FUNCTION
// #define KOKKOS_FUNCTION
// #endif

namespace kokkostbx {

template <typename NumType>
struct vector2 : public vector_base<vector2<NumType>, NumType, 2> {
    using vector_base = kokkostbx::vector_base<vector2<NumType>, NumType, 2>;

    vector2() = default;
    KOKKOS_FUNCTION vector2(NumType val) : vector_base(val){};
    KOKKOS_FUNCTION vector2(NumType arr[]) : vector_base(arr){};
    KOKKOS_FUNCTION vector2(const vector_base& vec) : vector_base(vec){};

    KOKKOS_FUNCTION vector2(NumType x, NumType y) : vector_base() {
        vector_base::data[0] = x;
        vector_base::data[1] = y;
    }

    // decided against using properties, as this would increase the size of the class
    KOKKOS_FUNCTION NumType& x_val() { return vector_base::data[0]; }
    KOKKOS_FUNCTION NumType& y_val() { return vector_base::data[1]; }

    KOKKOS_FUNCTION NumType x_val() const { return vector_base::data[0]; }
    KOKKOS_FUNCTION NumType y_val() const { return vector_base::data[1]; }
};

}  // namespace kokkostbx

#endif
