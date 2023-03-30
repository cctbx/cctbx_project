#ifndef KOKKOSTBX_MATRIX_H
#define KOKKOSTBX_MATRIX_H

#include "kokkos_vector.h"

// #ifndef KOKKOS_FUNCTION
// #define KOKKOS_FUNCTION
// #endif

namespace kokkostbx {

template <typename Derived, typename NumType, size_t rank>
struct matrix_base : public vector_base<matrix_base<Derived, NumType, rank>, NumType, rank * rank> {
    static constexpr size_t size = rank * rank;

    // data stored in row-major format!
    using vector_base =
        kokkostbx::vector_base<matrix_base<Derived, NumType, rank>, NumType, rank * rank>;
    using vector = kokkostbx::vector<NumType, rank>;

    matrix_base() = default;
    KOKKOS_FUNCTION matrix_base(NumType val) : vector_base(val){};
    KOKKOS_FUNCTION matrix_base(NumType arr[]) : vector_base(arr){};

    static Derived KOKKOS_FUNCTION diagonal_matrix(NumType arr[]) {
        Derived diagonal;
        for (size_t i = 0; i < rank; ++i) {
            diagonal.data[i * (rank + 1)] = arr[i];
        }
        return diagonal;
    };

    KOKKOS_FUNCTION NumType& operator()(size_t row, size_t col) {
        assert(row < rank);
        assert(col < rank);
        return vector_base::data[row * rank + col];
    }

    KOKKOS_FUNCTION NumType operator()(size_t row, size_t col) const {
        assert(row < rank);
        assert(col < rank);
        return vector_base::data[row * rank + col];
    }

    template <typename VectorType>
    KOKKOS_FUNCTION VectorType dot(const VectorType v) const {
        static_assert(v.get_size() == get_rank(), "Vector length must be equal to matrix rank!");
        VectorType result;
        for (size_t i = 0; i < rank; ++i) {
            for (size_t j = 0; j < rank; ++j) {
                result[i] += vector_base::data[i * rank + j] * v.data[j];
            }
        }
        return result;
    }

    KOKKOS_FUNCTION Derived dot(const Derived v) const {
        Derived result;
        for (size_t i = 0; i < rank; ++i) {
            for (size_t j = 0; j < rank; ++j) {
                for (size_t k = 0; k < rank; ++k) {
                    result[i * rank + j] += vector_base::data[i * rank + k] * v.data[k * rank + j];
                }
            }
        }
        return result;
    }

    static constexpr KOKKOS_FUNCTION size_t get_rank() { return rank; }

    KOKKOS_FUNCTION NumType trace() const {
        NumType result = 0;
        for (size_t i = 0; i < rank; ++i) {
            result += vector_base::data[i * rank + i];
        }
        return result;
    }

    KOKKOS_FUNCTION void print(const char name[]) const {
        printf("%s: ", name);
        for (size_t i = 0; i < size; ++i) {
            if (i % rank == 0) {
                printf("\n");
            }
            print_num(vector_base::data[i]);
        }
        printf("\n");
    }

    // KOKKOS_FUNCTION void transpose() {
    //     for (size_t i=0; i<(rank-1); ++i) {
    //         for (size_t j=(i+1); j<rank; ++j) {
    //             NumType temp = vector_base::data[i*rank+j];
    //             vector_base::data[i*rank+j] = vector_base::data[j*rank+i];
    //             vector_base::data[j*rank+i] = temp;
    //         }
    //     }
    // }

    KOKKOS_FUNCTION Derived transpose() const {
        Derived new_mat;
        for (size_t i = 0; i < rank; ++i) {
            for (size_t j = 0; j < rank; ++j) {
                new_mat[j * rank + i] = vector_base::data[i * rank + j];
            }
        }
        return new_mat;
    }
};

template <typename NumType, size_t rank>
struct matrix : public matrix_base<matrix<NumType, rank>, NumType, rank> {
    using matrix_base = kokkostbx::matrix_base<matrix<NumType, rank>, NumType, rank>;

    matrix() = default;
    KOKKOS_FUNCTION matrix(NumType val) : matrix_base(val){};
    KOKKOS_FUNCTION matrix(NumType arr[]) : matrix_base(arr){};
};

}  // namespace kokkostbx

#endif
