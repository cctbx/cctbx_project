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
    using VectorType = kokkostbx::vector<NumType, rank>;

    matrix_base() = default;
    KOKKOS_FUNCTION matrix_base(NumType val) : vector_base(val) {};
    KOKKOS_FUNCTION matrix_base(NumType arr[]) : vector_base(arr) {};

    static Derived KOKKOS_FUNCTION diagonal_matrix(NumType arr[]) {
        Derived diagonal;
        for (size_t i = 0; i < rank; ++i) {
            diagonal.data[i * (rank + 1)] = arr[i];
        }
        return diagonal;
    };

    // access
    KOKKOS_FUNCTION NumType& get(size_t row, size_t col) {
        assert(row < rank);
        assert(col < rank);
        return vector_base::data[row * rank + col];
    }

    KOKKOS_FUNCTION NumType get(size_t row, size_t col) const {
        assert(row < rank);
        assert(col < rank);
        return vector_base::data[row * rank + col];
    }

    KOKKOS_FUNCTION NumType& operator()(size_t row, size_t col) {
        return get(row, col);
    }

    KOKKOS_FUNCTION NumType operator()(size_t row, size_t col) const {
        return get(row, col);
    }

    // Override multiplication
    KOKKOS_FUNCTION friend Derived operator*(const Derived& lhs, const Derived& rhs) {
        Derived prod = lhs;
        prod *= rhs;
        return prod;
    }

    KOKKOS_FUNCTION friend VectorType operator*(const Derived& lhs, const VectorType& rhs) {
        return lhs.dot(rhs);
    }

    KOKKOS_FUNCTION friend Derived operator*(const Derived& lhs, const NumType& rhs) {
        Derived prod = lhs;
        prod *= rhs;
        return prod;
    }

    KOKKOS_FUNCTION friend Derived operator*(const NumType& lhs, const Derived& rhs) {
        return rhs * lhs;
    }

    KOKKOS_FUNCTION void operator*=(const Derived& rhs) {
        Derived prod = this->dot(rhs);
        for (size_t i = 0; i < size; ++i) {
            vector_base::data[i] = prod[i];
        }
    }

    KOKKOS_FUNCTION void operator*=(const NumType& rhs) {
        vector_base::operator*=(rhs);
    }

    KOKKOS_FUNCTION VectorType dot(const VectorType v) const {
        static_assert(v.get_size() == get_rank(), "Vector length must be equal to matrix rank!");
        VectorType result;
        for (size_t i = 0; i < rank; ++i) {
            for (size_t j = 0; j < rank; ++j) {
                result[i] += get(i, j) * v.data[j];
            }
        }
        return result;
    }

    KOKKOS_FUNCTION Derived dot(const Derived v) const {
        Derived result;
        for (size_t i = 0; i < rank; ++i) {
            for (size_t j = 0; j < rank; ++j) {
                for (size_t k = 0; k < rank; ++k) {
                    result(i, j) += get(i, k) * v(k, j);
                }
            }
        }
        return result;
    }

    static constexpr KOKKOS_FUNCTION size_t get_rank() { return rank; }

    KOKKOS_FUNCTION NumType trace() const {
        NumType result = 0;
        for (size_t i = 0; i < rank; ++i) {
            result += get(i, i);
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

    KOKKOS_FUNCTION Derived transpose() const {
        Derived new_mat;
        for (size_t i = 0; i < rank; ++i) {
            for (size_t j = 0; j < rank; ++j) {
                new_mat(j, i) = get(i, j);
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
