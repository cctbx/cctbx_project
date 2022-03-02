#ifndef SIMTBX_KOKKOS_VEC3_H
#define SIMTBX_KOKKOS_VEC3_H

#include <cmath>
#include <iostream>
#include <Kokkos_Core.hpp>

namespace simtbx { namespace kokkos {

    template <typename NumType>
    struct vec3 {
        NumType x = 0;
        NumType y = 0;
        NumType z = 0;

        // CONSTRUCTOR
        KOKKOS_FUNCTION vec3() = default;
        
        KOKKOS_FUNCTION vec3(NumType a) : x(a), y(a), z(a) {}

        KOKKOS_FUNCTION vec3(NumType a[3]) : x(a[0]), y(a[1]), z(a[2]) {}
        
        KOKKOS_FUNCTION vec3(NumType a, NumType b, NumType c) : x(a), y(b), z(c) {}

        // OPERATORS
        // streaming
        friend std::ostream& operator<< (std::ostream &os, const vec3<NumType>& v) {
            os << v.x << ", " << v.y << ", " << v.z;
            return os;
        }

        // access
        KOKKOS_INLINE_FUNCTION NumType& operator[](const int index) {
            switch(index) {
                case 0: return x;
                case 1: return y;
                case 2: return z;
                default:
                    ::Kokkos::abort("Can't access kokkos::vec3, index outside [0,1,2].");
            }
        }

        // addition
        KOKKOS_INLINE_FUNCTION friend vec3<NumType> operator+(const vec3<NumType>& lhs, const vec3<NumType>& rhs) {
            return vec3<NumType>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
        }
        KOKKOS_INLINE_FUNCTION friend vec3<NumType> operator+(const vec3<NumType>& lhs, NumType& rhs) {
            return vec3<NumType>(lhs.x+rhs, lhs.y+rhs, lhs.z+rhs);
        }
        KOKKOS_INLINE_FUNCTION friend vec3<NumType> operator+(NumType lhs, const vec3<NumType>& rhs) {
            return rhs + lhs;
        }
        
        KOKKOS_INLINE_FUNCTION void operator+=(const vec3<NumType>& v) {
            x += v.x;
            y += v.y;
            z += v.z;
        }

        KOKKOS_INLINE_FUNCTION void operator+=(const NumType& v) {
            x += v;
            y += v;
            z += v;
        }

        // subtraction
        KOKKOS_INLINE_FUNCTION vec3<NumType> operator-() const {
            return vec3<NumType>(-x, -y, -z);
        }

        KOKKOS_INLINE_FUNCTION vec3<NumType> operator-(const vec3<NumType>& v) const {
            return vec3<NumType>(x-v.x, y-v.y, z-v.z);
        }

        KOKKOS_INLINE_FUNCTION vec3<NumType> operator-(const NumType& v) const {
            return vec3<NumType>(x-v, y-v, z-v);
        }

        KOKKOS_INLINE_FUNCTION void operator-=(const vec3<NumType>& v) {
            x -= v.x;
            y -= v.y;
            z -= v.z;
        }

        KOKKOS_INLINE_FUNCTION void operator-=(const NumType& v) {
            x -= v;
            y -= v;
            z -= v;
        }

        // multiplication
        KOKKOS_INLINE_FUNCTION friend vec3<NumType> operator*(const vec3<NumType>& lhs, const vec3<NumType>& rhs) {
            return vec3<NumType>(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
        }
        KOKKOS_INLINE_FUNCTION friend vec3<NumType> operator*(const vec3<NumType>& lhs, NumType& rhs) {
            return vec3<NumType>(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
        }
        KOKKOS_INLINE_FUNCTION friend vec3<NumType> operator*(NumType lhs, const vec3<NumType>& rhs) {
            return rhs * lhs;
        }

        KOKKOS_INLINE_FUNCTION void operator*=(const vec3<NumType>& v) {
            x *= v.x;
            y *= v.y;
            z *= v.z;
        }

        KOKKOS_INLINE_FUNCTION void operator*=(const NumType& v) {
            x *= v;
            y *= v;
            z *= v;
        }

        // division
        KOKKOS_INLINE_FUNCTION vec3<NumType> operator/(const vec3<NumType>& v) const {
            return vec3<NumType>(x/v.x, y/v.y, z/v.z);
        }

        KOKKOS_INLINE_FUNCTION vec3<NumType> operator/(const NumType& v) const {
            return vec3<NumType>(x/v, y/v, z/v);
        }

        KOKKOS_INLINE_FUNCTION void operator/=(const vec3<NumType>& v) {
            x /= v.x;
            y /= v.y;
            z /= v.z;
        }

        KOKKOS_INLINE_FUNCTION void operator/=(const NumType& v) {
            x /= v;
            y /= v;
            z /= v;
        }

        // METHODS
        KOKKOS_INLINE_FUNCTION void zero() {
            x = 0;
            z = 0;
            y = 0;
        }

        KOKKOS_INLINE_FUNCTION NumType length_sqr() const {
            return x*x + y*y + z*z;
        }

        KOKKOS_INLINE_FUNCTION NumType length() const {
            return std::sqrt(length_sqr());
        }

        KOKKOS_INLINE_FUNCTION NumType dot(const vec3<NumType>& v) const {
            return x*v.x + y*v.y + z*v.z;
        }

        KOKKOS_INLINE_FUNCTION vec3<NumType> cross(const vec3<NumType>& v) const {
            return vec3<NumType>(y*v.z - z*v.y, 
                                z*v.x - x*v.z,
                                x*v.y - y*v.x);
        }

        KOKKOS_INLINE_FUNCTION void normalize() {
            NumType l = length();
            if (l>0) {
                x /= l;
                y /= l;
                z /= l;
            }
        }

        KOKKOS_INLINE_FUNCTION vec3<NumType> get_unit_vector() const {
            NumType l = length();
            if (l>0) {
                return vec3<NumType>(x/l, y/l, z/l);
            } else {
                return vec3<NumType>(0, 0, 0);
            }
        }

        // rotate a point about a unit vec3 axis
        KOKKOS_INLINE_FUNCTION vec3<NumType> rotate_around_axis(const vec3<NumType>& axis, NumType angle) const {
            NumType sinphi = std::sin(angle);
            NumType cosphi = std::cos(angle);
            NumType dot_factor = axis.dot(*this) * (1.0-cosphi);
            
            vec3<NumType> vector_rot = axis.cross(*this) * sinphi;
            vector_rot += axis * dot_factor;
            vector_rot += (*this) * cosphi;
            return vector_rot;
        }
        
    };
} } // namespace scitbx::kokkos




#endif // SIMTBX_KOKKOS_VEC3_H
