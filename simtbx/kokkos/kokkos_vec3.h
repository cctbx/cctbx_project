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

    // polarization factor
    /*template <typename NumType>
    KOKKOS_FUNCTION CUDAREAL polarization_factor(CUDAREAL kahn_factor, const vec3<NumType>& incident, const vec3<NumType>& diffracted, const vec3<NumType>& axis) {
        CUDAREAL cos2theta, cos2theta_sqr, sin2theta_sqr;
        CUDAREAL psi = 0.0;
        CUDAREAL E_in[4], B_in[4], E_out[4], B_out[4];

        //  these are already unitized before entering this loop. Optimize this out.
        //        unitize(incident, incident);
        //        unitize(diffracted, diffracted);

        // component of diffracted unit vector along incident beam unit vector
        cos2theta = dot_product(incident, diffracted);
        cos2theta_sqr = cos2theta * cos2theta;
        sin2theta_sqr = 1 - cos2theta_sqr;

        if (kahn_factor != 0.0) {
                // tricky bit here is deciding which direciton the E-vector lies in for each source
                // here we assume it is closest to the "axis" defined above

                CUDAREAL unitAxis[] = { axis(0), axis(1), axis(2), axis(3) };
                // this is already unitized. Optimize this out.
                unitize(unitAxis, unitAxis);

                // cross product to get "vertical" axis that is orthogonal to the cannonical "polarization"
                cross_product(unitAxis, incident, B_in);
                // make it a unit vector
                unitize(B_in, B_in);

                // cross product with incident beam to get E-vector direction
                cross_product(incident, B_in, E_in);
                // make it a unit vector
                unitize(E_in, E_in);

                // get components of diffracted ray projected onto the E-B plane
                E_out[0] = dot_product(diffracted, E_in);
                B_out[0] = dot_product(diffracted, B_in);

                // compute the angle of the diffracted ray projected onto the incident E-B plane
                psi = -atan2(B_out[0], E_out[0]);
        }

        // correction for polarized incident beam
        return 0.5 * (1.0 + cos2theta_sqr - kahn_factor * cos(2 * psi) * sin2theta_sqr);
    }*/
} } // namespace scitbx::kokkos




#endif // SIMTBX_KOKKOS_VEC3_H
