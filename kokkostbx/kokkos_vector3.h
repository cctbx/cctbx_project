#ifndef KOKKOSTBX_VECTOR3_H
#define KOKKOSTBX_VECTOR3_H

#include "kokkos_vector.h"

// Bandaid since parts of diffBragg don't support Kokkos yet
// #ifndef KOKKOS_FUNCTION
// #define KOKKOS_FUNCTION
// #endif

// #ifdef KOKKOS_CORE_HPP
//     template <typename T> KOKKOS_FUNCTION T sin_func(T x) { return
//     ::Kokkos::sin(x); } template <typename T> KOKKOS_FUNCTION T cos_func(T x) {
//     return ::Kokkos::cos(x); }
// #else
//     #include <cmath>
//     template <typename T> KOKKOS_FUNCTION T sin_func(T x) { return sin(x); }
//     template <typename T> KOKKOS_FUNCTION T cos_func(T x) { return cos(x); }
// #endif

namespace kokkostbx {

template <typename NumType>
//struct vector3 : public vector_base<vector3<NumType>, NumType, 3> {
struct vector3 : public vector<NumType, 3> {
    //using vector_base = kokkostbx::vector_base<vector<NumType>, NumType, 3>;
    using vector_base = kokkostbx::vector<NumType, 3>;

    vector3() = default;
    KOKKOS_INLINE_FUNCTION vector3(NumType val) : vector_base(val){};
    KOKKOS_INLINE_FUNCTION vector3(NumType arr[]) : vector_base(arr){};
    KOKKOS_INLINE_FUNCTION vector3(const vector_base& vec) : vector_base(vec){};

    KOKKOS_INLINE_FUNCTION vector3(NumType x, NumType y, NumType z) : vector_base() {
        vector_base::data[0] = x;
        vector_base::data[1] = y;
        vector_base::data[2] = z;
    }

    // decided against using properties, as this would increase the size of the class
    KOKKOS_INLINE_FUNCTION NumType& x_val() { return vector_base::data[0]; }
    KOKKOS_INLINE_FUNCTION NumType& y_val() { return vector_base::data[1]; }
    KOKKOS_INLINE_FUNCTION NumType& z_val() { return vector_base::data[2]; }

    KOKKOS_INLINE_FUNCTION NumType x_val() const { return vector_base::data[0]; }
    KOKKOS_INLINE_FUNCTION NumType y_val() const { return vector_base::data[1]; }
    KOKKOS_INLINE_FUNCTION NumType z_val() const { return vector_base::data[2]; }

    KOKKOS_INLINE_FUNCTION vector3<NumType> cross(const vector3<NumType>& v) const {
        vector3<NumType> cross_vector{};
        cross_vector.x_val() = y_val() * v.z_val() - z_val() * v.y_val();
        cross_vector.y_val() = z_val() * v.x_val() - x_val() * v.z_val();
        cross_vector.z_val() = x_val() * v.y_val() - y_val() * v.x_val();

        return cross_vector;
    }

    // rotate a point around a unit vector3 axis
    KOKKOS_INLINE_FUNCTION vector3<NumType> rotate_around_axis(const vector3<NumType>& axis, NumType angle)
        const {
        // NumType sinphi = sin_func(angle);
        // NumType cosphi = cos_func(angle);
        const NumType sinphi = ::Kokkos::sin(angle);
        const NumType cosphi = ::Kokkos::cos(angle);

        const NumType dot_factor = axis.dot(*this) * (1.0 - cosphi);

        vector3<NumType> vector_rot = axis.cross(*this) * sinphi;
        vector_rot += axis * dot_factor;
        vector_rot += (*this) * cosphi;
        return vector_rot;
    }
};

}  // namespace kokkostbx

#endif
