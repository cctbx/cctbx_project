#ifndef KOKKOSTBX_MATRIX_H
#define KOKKOSTBX_MATRIX_H

#include <Kokkos_Core.hpp>

#include "kokkos_vector.h"

namespace kokkostbx {

    template <typename Derived, typename NumType, size_t rank>
    struct matrix_base : public vector_base<matrix_base<Derived, NumType, rank>, NumType, rank*rank> { 

        const size_t size = rank*rank;

        // data stored in row-major format!
        using vector_base = kokkostbx::vector_base<matrix_base<Derived, NumType, rank>, NumType, rank*rank>;
        using vector = kokkostbx::vector<NumType, rank>;

        KOKKOS_FUNCTION matrix_base() = default;
        KOKKOS_FUNCTION matrix_base(NumType val) : vector_base(val) { };
        KOKKOS_FUNCTION matrix_base(NumType arr[]) : vector_base(arr) { };

        static Derived KOKKOS_FUNCTION diagonal_matrix(NumType arr[]) {
            Derived diagonal; 
            for (size_t i=0; i<rank; ++i) {
                diagonal.data[i*(rank+1)] = arr[i];
            }
            return diagonal;
        };

        template <typename VectorType>
        KOKKOS_FUNCTION VectorType dot(const VectorType v) const {
            static_assert(v.get_size() == get_rank(), "Vector length must be equal to matrix rank!");
            VectorType result;
            for (size_t i=0; i<rank; ++i) {
                for (size_t j=0; j<rank; ++j) {
                    result[i] += vector_base::data[i*rank+j] * v.data[j];
                }
            }
            return result;
        }

        static constexpr KOKKOS_FUNCTION size_t get_rank() {
            return rank;
        }

        KOKKOS_FUNCTION void print(const char name[]) const {
            printf("%s: ", name);
            for (size_t i=0; i<size; ++i) {
                if (i%rank==0) { printf("\n"); }
                print_num(vector_base::data[i]);
            }
            printf("\n");
        } 
    };

    template <typename NumType, size_t rank>
    struct matrix : public matrix_base<matrix<NumType, rank>, NumType, rank> { 

        using matrix_base = kokkostbx::matrix_base<matrix<NumType, rank>, NumType, rank>;

        matrix() = default;
        KOKKOS_FUNCTION matrix(NumType val) : matrix_base(val) { };
        KOKKOS_FUNCTION matrix(NumType arr[]) : matrix_base(arr) { };
    };

}

#endif
