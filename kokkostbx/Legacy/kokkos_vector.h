#ifndef KOKKOSTBX_VECTOR_H
#define KOKKOSTBX_VECTOR_H

#include <iostream>

namespace kokkostbx {

    template <typename Derived, typename NumType, size_t size>
    struct vector_base {

        NumType data[size] = {};

        // CONSTRUCTOR
        KOKKOS_FUNCTION vector_base() = default;

        KOKKOS_FUNCTION vector_base(NumType val) {
            for (NumType& d : data) { d = val; }
        }

        KOKKOS_FUNCTION vector_base(NumType arr[]) {
            for (int i=0; i<size; ++i) {
                data[i] = arr[i];
            }
        }

        // OPERATORS
        // streaming
        friend std::ostream& operator<< (std::ostream &os, const vector_base<Derived, NumType, size>& v) {
            for (int i=0; i<size: ++i ) {
                if (i>0) { os << ", "; }
                os << data[i];
            }
            return os;
        }

        // access
        KOKKOS_INLINE_FUNCTION NumType& operator[](const int index) {
            return data[index];
        }

        // addition
        KOKKOS_INLINE_FUNCTION friend vector_base<Derived, NumType, size> operator+(
            const vector_base<Derived, NumType, size>& lhs,
            const vector_base<Derived, NumType, size>& rhs) {

            vector_base<Derived, NumType, size> sum = lhs;
            sum += rhs;
            return sum;
        }
        KOKKOS_INLINE_FUNCTION friend vector_base<Derived, NumType, size> operator+(
            const vector_base<Derived, NumType, size>& lhs,
            NumType& rhs) {

            vector_base<Derived, NumType, size> sum = lhs;
            sum += rhs;
            return sum;
        }
        KOKKOS_INLINE_FUNCTION friend vector_base<Derived, NumType, size> operator+(
            NumType lhs,
            const vector_base<Derived, NumType, size>& rhs) {
            return rhs + lhs;
        }

        KOKKOS_INLINE_FUNCTION void operator+=(const vector_base<Derived, NumType, size>& v) {
            for (int i=0; i<size; ++i) {
                data[i] += v[i];
            }
        }

        KOKKOS_INLINE_FUNCTION void operator+=(const NumType& v) {
            for (int i=0; i<size; ++i) {
                data[i] += v;
            }
        }
    }

    };



    template <typename NumType>
    struct vector3 {
        NumType x = 0;
        NumType y = 0;
        NumType z = 0;

        // CONSTRUCTOR
        KOKKOS_FUNCTION vector3() = default;

        KOKKOS_FUNCTION vector3(NumType a) : x(a), y(a), z(a) {}

        KOKKOS_FUNCTION vector3(NumType a[3]) : x(a[0]), y(a[1]), z(a[2]) {}

        KOKKOS_FUNCTION vector3(NumType a, NumType b, NumType c) : x(a), y(b), z(c) {}

        // OPERATORS
        // streaming
        friend std::ostream& operator<< (std::ostream &os, const vector3<NumType>& v) {
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
                    ::Kokkos::abort("Can't access kokkos::vector3, index outside [0,1,2].");
            }
        }

        // addition
        KOKKOS_INLINE_FUNCTION friend vector3<NumType> operator+(const vector3<NumType>& lhs, const vector3<NumType>& rhs) {
            return vector3<NumType>(lhs.x+rhs.x, lhs.y+rhs.y, lhs.z+rhs.z);
        }
        KOKKOS_INLINE_FUNCTION friend vector3<NumType> operator+(const vector3<NumType>& lhs, NumType& rhs) {
            return vector3<NumType>(lhs.x+rhs, lhs.y+rhs, lhs.z+rhs);
        }
        KOKKOS_INLINE_FUNCTION friend vector3<NumType> operator+(NumType lhs, const vector3<NumType>& rhs) {
            return rhs + lhs;
        }

        KOKKOS_INLINE_FUNCTION void operator+=(const vector3<NumType>& v) {
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
        KOKKOS_INLINE_FUNCTION vector3<NumType> operator-() const {
            return vector3<NumType>(-x, -y, -z);
        }

        KOKKOS_INLINE_FUNCTION vector3<NumType> operator-(const vector3<NumType>& v) const {
            return vector3<NumType>(x-v.x, y-v.y, z-v.z);
        }

        KOKKOS_INLINE_FUNCTION vector3<NumType> operator-(const NumType& v) const {
            return vector3<NumType>(x-v, y-v, z-v);
        }

        KOKKOS_INLINE_FUNCTION void operator-=(const vector3<NumType>& v) {
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
        KOKKOS_INLINE_FUNCTION friend vector3<NumType> operator*(const vector3<NumType>& lhs, const vector3<NumType>& rhs) {
            return vector3<NumType>(lhs.x*rhs.x, lhs.y*rhs.y, lhs.z*rhs.z);
        }
        KOKKOS_INLINE_FUNCTION friend vector3<NumType> operator*(const vector3<NumType>& lhs, NumType& rhs) {
            return vector3<NumType>(lhs.x*rhs, lhs.y*rhs, lhs.z*rhs);
        }
        KOKKOS_INLINE_FUNCTION friend vector3<NumType> operator*(NumType lhs, const vector3<NumType>& rhs) {
            return rhs * lhs;
        }

        KOKKOS_INLINE_FUNCTION void operator*=(const vector3<NumType>& v) {
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
        KOKKOS_INLINE_FUNCTION vector3<NumType> operator/(const vector3<NumType>& v) const {
            return vector3<NumType>(x/v.x, y/v.y, z/v.z);
        }

        KOKKOS_INLINE_FUNCTION vector3<NumType> operator/(const NumType& v) const {
            return vector3<NumType>(x/v, y/v, z/v);
        }

        KOKKOS_INLINE_FUNCTION void operator/=(const vector3<NumType>& v) {
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
            return ::Kokkos::sqrt(length_sqr());
        }

        KOKKOS_INLINE_FUNCTION NumType dot(const vector3<NumType>& v) const {
            return x*v.x + y*v.y + z*v.z;
        }

        KOKKOS_INLINE_FUNCTION vector3<NumType> cross(const vector3<NumType>& v) const {
            return vector3<NumType>(y*v.z - z*v.y,
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

        KOKKOS_INLINE_FUNCTION vector3<NumType> get_unit_vector() const {
            NumType l = length();
            if (l>0) {
                return vector3<NumType>(x/l, y/l, z/l);
            } else {
                return vector3<NumType>(0, 0, 0);
            }
        }

        // rotate a point about a unit vector3 axis
        KOKKOS_INLINE_FUNCTION vector3<NumType> rotate_around_axis(const vector3<NumType>& axis, NumType angle) const {
            NumType sinphi = ::Kokkos::sin(angle);
            NumType cosphi = ::Kokkos::cos(angle);
            NumType dot_factor = axis.dot(*this) * (1.0-cosphi);

            vector3<NumType> vector_rot = axis.cross(*this) * sinphi;
            vector_rot += axis * dot_factor;
            vector_rot += (*this) * cosphi;
            return vector_rot;
        }

        // rotate a vector using a 9-element unitary matrix
        KOKKOS_INLINE_FUNCTION vector3<NumType> rotate_matrix(const NumType * __restrict__ umat) {
            // for convenience, assign matrix x-y coordinate
            vector3<NumType> ux = {umat[0], umat[1], umat[2]};
            vector3<NumType> uy = {umat[3], umat[4], umat[5]};
            vector3<NumType> uz = {umat[6], umat[7], umat[8]};

            // rotate the vector
            NumType newx = ux.dot(*this);
            NumType newy = uy.dot(*this);
            NumType newz = uz.dot(*this);

            return vector3<NumType>(newx, newy, newz);
        }
    };

}

#endif




#endif
