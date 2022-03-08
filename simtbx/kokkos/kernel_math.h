#ifndef SIMTBX_KOKKOS_KERNEL_MATH_H
#define SIMTBX_KOKKOS_KERNEL_MATH_H

#ifndef CUDAREAL
#define CUDAREAL float
#endif

// cubic spline interpolation functions
KOKKOS_FUNCTION static CUDAREAL polint(const CUDAREAL *xa, const CUDAREAL *ya, CUDAREAL x);

// cubic spline interpolation functions
KOKKOS_FUNCTION CUDAREAL polint(const CUDAREAL *xa, const CUDAREAL *ya, CUDAREAL x) {
    auto interpolate = [x](CUDAREAL yt, CUDAREAL xt, CUDAREAL a, CUDAREAL b, CUDAREAL c) {
        return (x - a) * (x - b) * (x - c) * yt / ((xt - a) * (xt - b) * (xt - c));
    };
    CUDAREAL t0 = interpolate(ya[0], xa[0], xa[1], xa[2], xa[3]);
    CUDAREAL t1 = interpolate(ya[1], xa[1], xa[0], xa[2], xa[3]);
    CUDAREAL t2 = interpolate(ya[2], xa[2], xa[0], xa[1], xa[3]);
    CUDAREAL t3 = interpolate(ya[3], xa[3], xa[0], xa[1], xa[2]);
    return t0 + t1 + t2 + t3;
}

namespace simtbx { namespace kokkos { 
    // polarization factor
    template <typename NumType>
    KOKKOS_FUNCTION CUDAREAL polarization_factor(NumType kahn_factor, const vec3<NumType>& incident, const vec3<NumType>& diffracted, const vec3<NumType>& axis) {
        NumType psi = 0.0;

        // component of diffracted unit vector along incident beam unit vector
        NumType cos2theta = incident.dot(diffracted);
        NumType cos2theta_sqr = cos2theta * cos2theta;
        NumType sin2theta_sqr = 1 - cos2theta_sqr;

        if (kahn_factor != 0.0) {
            // tricky bit here is deciding which direction the E-vector lies in for each source
            // here we assume it is closest to the "axis" defined above

            // cross product to get "vertical" axis that is orthogonal to the canonical "polarization"
            vec3<NumType> unitAxis = axis.get_unit_vector();
            vec3<NumType> B_in = unitAxis.cross(incident);
            B_in.normalize();

            // cross product with incident beam to get E-vector direction
            vec3<NumType> E_in = incident.cross(B_in);
            E_in.normalize();

            // get components of diffracted ray projected onto the E-B plane
            CUDAREAL E_out = diffracted.dot(E_in);
            CUDAREAL B_out = diffracted.dot(B_in);

            // compute the angle of the diffracted ray projected onto the incident E-B plane
            psi = -atan2(B_out, E_out);
        }

        // correction for polarized incident beam
        return 0.5 * (1.0 + cos2theta_sqr - kahn_factor * cos(2 * psi) * sin2theta_sqr);
    }
} } // namespace simtbx::kokkos

// Fourier transform of a grating
KOKKOS_INLINE_FUNCTION CUDAREAL sincg(CUDAREAL x, CUDAREAL N) {
    if (x != 0.0) {
        return sin(x * N) / sin(x);
    }
    return N;
}

#ifndef sinpi
#define sinpi(x) sin(x * M_PI)
#endif

KOKKOS_INLINE_FUNCTION CUDAREAL sincgrad(CUDAREAL x, CUDAREAL N) {
    if (x != 0.0) {
        return sinpi(x * N) / sinpi(x);
    }
    return N;
}

// Fourier transform of a sphere
KOKKOS_INLINE_FUNCTION CUDAREAL sinc3(CUDAREAL x) {
    if (x != 0.0) {
        return 3.0 * (sin(x) / x - cos(x)) / (x * x);
    }
    return 1.0;
}

#endif // SIMTBX_KOKKOS_KERNEL_MATH_H
