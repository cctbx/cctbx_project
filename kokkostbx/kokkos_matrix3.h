#ifndef KOKKOSTBX_MATRIX3_H
#define KOKKOSTBX_MATRIX3_H

#include "kokkos_matrix.h"

// #ifndef KOKKOS_FUNCTION
// #define KOKKOS_FUNCTION
// #endif

namespace kokkostbx {

template <typename NumType>
struct matrix3 : public matrix_base<matrix3<NumType>, NumType, 3> {
    using matrix_base = kokkostbx::matrix_base<matrix3<NumType>, NumType, 3>;

    matrix3() = default;
    KOKKOS_FUNCTION matrix3(NumType val) : matrix_base(val){};
    KOKKOS_FUNCTION matrix3(NumType arr[]) : matrix_base(arr){};

    KOKKOS_FUNCTION matrix3(NumType x, NumType y, NumType z) : matrix_base() {
        matrix_base::data[0 * 3 + 0] = x;
        matrix_base::data[1 * 3 + 1] = y;
        matrix_base::data[2 * 3 + 2] = z;
    }

    KOKKOS_FUNCTION Derived determinant() const {
        NumType result = 0;

        result += vector_base[0 * 3 + 0] * vector_base[1 * 3 + 1] * vector_base[2 * 3 + 2];
        result -= vector_base[0 * 3 + 2] * vector_base[1 * 3 + 1] * vector_base[2 * 3 + 0];
        result += vector_base[0 * 3 + 1] * vector_base[1 * 3 + 2] * vector_base[2 * 3 + 0];
        result -= vector_base[0 * 3 + 0] * vector_base[1 * 3 + 2] * vector_base[2 * 3 + 1];
        result += vector_base[0 * 3 + 2] * vector_base[1 * 3 + 0] * vector_base[2 * 3 + 1];        
        result -= vector_base[0 * 3 + 1] * vector_base[1 * 3 + 0] * vector_base[2 * 3 + 2];
        return result;
    }
};
}  // namespace kokkostbx

#endif
