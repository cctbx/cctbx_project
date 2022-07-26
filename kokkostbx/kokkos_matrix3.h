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

    KOKKOS_FUNCTION matrix3(
        NumType a,
        NumType b,
        NumType c,
        NumType d,
        NumType e,
        NumType f,
        NumType g,
        NumType h,
        NumType i)
        : matrix_base() {
        matrix_base::data[0 * 3 + 0] = a;
        matrix_base::data[0 * 3 + 1] = b;
        matrix_base::data[0 * 3 + 2] = c;
        matrix_base::data[1 * 3 + 0] = d;
        matrix_base::data[1 * 3 + 1] = e;
        matrix_base::data[1 * 3 + 2] = f;
        matrix_base::data[2 * 3 + 0] = g;
        matrix_base::data[2 * 3 + 1] = h;
        matrix_base::data[2 * 3 + 2] = i;
    }

    KOKKOS_FUNCTION matrix3 inverse() const {
        matrix3 result;
        // Full matrix:
        //  a b c
        //  d e f
        //  g h i
        const NumType a = matrix_base::data[0 * 3 + 0];
        const NumType b = matrix_base::data[0 * 3 + 1];
        const NumType c = matrix_base::data[0 * 3 + 2];
        const NumType d = matrix_base::data[1 * 3 + 0];
        const NumType e = matrix_base::data[1 * 3 + 1];
        const NumType f = matrix_base::data[1 * 3 + 2];
        const NumType g = matrix_base::data[2 * 3 + 0];
        const NumType h = matrix_base::data[2 * 3 + 1];
        const NumType i = matrix_base::data[2 * 3 + 2];

        result[0 * 3 + 0] = e * i - f * h;
        result[0 * 3 + 1] = c * h - b * i;
        result[0 * 3 + 2] = b * f - c * e;
        result[1 * 3 + 0] = f * g - d * i;
        result[1 * 3 + 1] = a * i - c * g;
        result[1 * 3 + 2] = c * d - a * f;
        result[2 * 3 + 0] = d * h - e * g;
        result[2 * 3 + 1] = b * g - a * h;
        result[2 * 3 + 2] = a * e - b * d;

        result /= determinant();
        return result;
    }

    KOKKOS_FUNCTION NumType determinant() const {
        NumType result = 0;

        result += matrix_base::data[0 * 3 + 0] * matrix_base::data[1 * 3 + 1] *
                  matrix_base::data[2 * 3 + 2];
        result -= matrix_base::data[0 * 3 + 2] * matrix_base::data[1 * 3 + 1] *
                  matrix_base::data[2 * 3 + 0];
        result += matrix_base::data[0 * 3 + 1] * matrix_base::data[1 * 3 + 2] *
                  matrix_base::data[2 * 3 + 0];
        result -= matrix_base::data[0 * 3 + 0] * matrix_base::data[1 * 3 + 2] *
                  matrix_base::data[2 * 3 + 1];
        result += matrix_base::data[0 * 3 + 2] * matrix_base::data[1 * 3 + 0] *
                  matrix_base::data[2 * 3 + 1];
        result -= matrix_base::data[0 * 3 + 1] * matrix_base::data[1 * 3 + 0] *
                  matrix_base::data[2 * 3 + 2];
        return result;
    }
};
}  // namespace kokkostbx

#endif
