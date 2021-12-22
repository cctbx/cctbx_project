#ifndef SIMTBX_KOKKOS_KERNEL_MATH_H
#define SIMTBX_KOKKOS_KERNEL_MATH_H

#ifndef CUDAREAL
#define CUDAREAL float
#endif

// cubic spline interpolation functions
KOKKOS_FUNCTION static void polint(const CUDAREAL *xa, const CUDAREAL *ya, CUDAREAL x, CUDAREAL *y);
// vector inner product where vector magnitude is 0th element
KOKKOS_FUNCTION static CUDAREAL dot_product(const CUDAREAL *x, const CUDAREAL *y);
// KOKKOS_INLINE_FUNCTION static CUDAREAL dot_product_ldg(const CUDAREAL * __restrict__ x, CUDAREAL * y);
// make a unit vector pointing in same direction and report magnitude (both args can be same vector)
KOKKOS_FUNCTION static CUDAREAL unitize(CUDAREAL * vector, CUDAREAL *new_unit_vector);
// polarization factor from vectors
KOKKOS_FUNCTION static CUDAREAL polarization_factor(CUDAREAL kahn_factor, CUDAREAL *incident, CUDAREAL *diffracted, const vector_cudareal_t axis);
// vector cross product where vector magnitude is 0th element
KOKKOS_INLINE_FUNCTION static void cross_product(CUDAREAL *x, CUDAREAL *y, CUDAREAL *z);
/* rotate a 3-vector about a unit vector axis */
KOKKOS_INLINE_FUNCTION static CUDAREAL *rotate_axis(const vector_cudareal_t v, CUDAREAL * newv, const vector_cudareal_t axis, const CUDAREAL phi);
KOKKOS_INLINE_FUNCTION static CUDAREAL *rotate_axis(const CUDAREAL * v, CUDAREAL * newv, const vector_cudareal_t axis, const CUDAREAL phi);
/* rotate a 3-vector using a 9-element unitary matrix */
KOKKOS_INLINE_FUNCTION static void rotate_umat(CUDAREAL * v, CUDAREAL *newv, const CUDAREAL * __restrict__ umat);
// measure magnitude of vector and put it in 0th element
KOKKOS_FUNCTION static void magnitude(CUDAREAL *vector);

// cubic spline interpolation functions
KOKKOS_FUNCTION void polint(const CUDAREAL *xa, const CUDAREAL *ya, CUDAREAL x, CUDAREAL *y) {
        CUDAREAL t0 = (x - xa[1]) * (x - xa[2]) * (x - xa[3]) * ya[0] / ((xa[0] - xa[1]) * (xa[0] - xa[2]) * (xa[0] - xa[3]));
        CUDAREAL t1 = (x - xa[0]) * (x - xa[2]) * (x - xa[3]) * ya[1] / ((xa[1] - xa[0]) * (xa[1] - xa[2]) * (xa[1] - xa[3]));
        CUDAREAL t2 = (x - xa[0]) * (x - xa[1]) * (x - xa[3]) * ya[2] / ((xa[2] - xa[0]) * (xa[2] - xa[1]) * (xa[2] - xa[3]));
        CUDAREAL t3 = (x - xa[0]) * (x - xa[1]) * (x - xa[2]) * ya[3] / ((xa[3] - xa[0]) * (xa[3] - xa[1]) * (xa[3] - xa[2]));
        *y = t0 + t1 + t2 + t3;
}

// vector inner product where vector magnitude is 0th element
KOKKOS_FUNCTION CUDAREAL dot_product(const CUDAREAL *x, const CUDAREAL *y) {
        return x[1] * y[1] + x[2] * y[2] + x[3] * y[3];
}
// KOKKOS_INLINE_FUNCTION CUDAREAL dot_product_ldg(const CUDAREAL * __restrict__ x, CUDAREAL * y) {
//         return __ldg(&x[1]) * y[1] + __ldg(&x[2]) * y[2] + __ldg(&x[3]) * y[3];
// }

// make provided vector a unit vector
KOKKOS_FUNCTION CUDAREAL unitize(CUDAREAL * vector, CUDAREAL * new_unit_vector) {

        CUDAREAL v1 = vector[1];
        CUDAREAL v2 = vector[2];
        CUDAREAL v3 = vector[3];

#ifdef __CUDA_ARCH__
        CUDAREAL mag = norm3d(v1, v2, v3);
#else
        CUDAREAL mag = sqrt(v1 * v1 + v2 * v2 + v3 * v3);
#endif

        if (mag != 0.0) {
                new_unit_vector[0] = mag;
                new_unit_vector[1] = v1 / mag;
                new_unit_vector[2] = v2 / mag;
                new_unit_vector[3] = v3 / mag;
        } else {
                new_unit_vector[0] = 0.0;
                new_unit_vector[1] = 0.0;
                new_unit_vector[2] = 0.0;
                new_unit_vector[3] = 0.0;
        }
        return mag;
}

// polarization factor
KOKKOS_FUNCTION CUDAREAL polarization_factor(CUDAREAL kahn_factor, CUDAREAL *incident, CUDAREAL *diffracted, const vector_cudareal_t axis) {
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
}

// vector cross product where vector magnitude is 0th element
KOKKOS_INLINE_FUNCTION void cross_product(CUDAREAL *x, CUDAREAL *y, CUDAREAL *z) {
        z[1] = x[2] * y[3] - x[3] * y[2];
        z[2] = x[3] * y[1] - x[1] * y[3];
        z[3] = x[1] * y[2] - x[2] * y[1];
        z[0] = 0.0;
}

/* rotate a point about a unit vector axis */
KOKKOS_INLINE_FUNCTION CUDAREAL *rotate_axis(const vector_cudareal_t v, CUDAREAL * newv, const vector_cudareal_t axis, const CUDAREAL phi) {

        const CUDAREAL sinphi = sin(phi);
        const CUDAREAL cosphi = cos(phi);
        const CUDAREAL a1 = axis(1);
        const CUDAREAL a2 = axis(2);
        const CUDAREAL a3 = axis(3);
        const CUDAREAL v1 = v(1);
        const CUDAREAL v2 = v(2);
        const CUDAREAL v3 = v(3);
        const CUDAREAL dot = (a1 * v1 + a2 * v2 + a3 * v3) * (1.0 - cosphi);

        newv[1] = a1 * dot + v1 * cosphi + (-a3 * v2 + a2 * v3) * sinphi;
        newv[2] = a2 * dot + v2 * cosphi + (+a3 * v1 - a1 * v3) * sinphi;
        newv[3] = a3 * dot + v3 * cosphi + (-a2 * v1 + a1 * v2) * sinphi;

        return newv;
}

/* rotate a point about a unit vector axis */
KOKKOS_INLINE_FUNCTION CUDAREAL *rotate_axis(const CUDAREAL * v, CUDAREAL * newv, const vector_cudareal_t axis, const CUDAREAL phi) {

        const CUDAREAL sinphi = sin(phi);
        const CUDAREAL cosphi = cos(phi);
        const CUDAREAL a1 = axis(1);
        const CUDAREAL a2 = axis(2);
        const CUDAREAL a3 = axis(3);
        const CUDAREAL v1 = v[1];
        const CUDAREAL v2 = v[2];
        const CUDAREAL v3 = v[3];
        const CUDAREAL dot = (a1 * v1 + a2 * v2 + a3 * v3) * (1.0 - cosphi);

        newv[1] = a1 * dot + v1 * cosphi + (-a3 * v2 + a2 * v3) * sinphi;
        newv[2] = a2 * dot + v2 * cosphi + (+a3 * v1 - a1 * v3) * sinphi;
        newv[3] = a3 * dot + v3 * cosphi + (-a2 * v1 + a1 * v2) * sinphi;

        return newv;
}

/* rotate a vector using a 9-element unitary matrix */
KOKKOS_INLINE_FUNCTION void rotate_umat(CUDAREAL * v, CUDAREAL *newv, const CUDAREAL * __restrict__ umat) {

        /* for convenience, assign matrix x-y coordinate */
        CUDAREAL uxx = umat[0];
        CUDAREAL uxy = umat[1];
        CUDAREAL uxz = umat[2];
        CUDAREAL uyx = umat[3];
        CUDAREAL uyy = umat[4];
        CUDAREAL uyz = umat[5];
        CUDAREAL uzx = umat[6];
        CUDAREAL uzy = umat[7];
        CUDAREAL uzz = umat[8];
        CUDAREAL v1 = v[1];
        CUDAREAL v2 = v[2];
        CUDAREAL v3 = v[3];

        /* rotate the vector (x=1,y=2,z=3) */
        newv[1] = uxx * v1 + uxy * v2 + uxz * v3;
        newv[2] = uyx * v1 + uyy * v2 + uyz * v3;
        newv[3] = uzx * v1 + uzy * v2 + uzz * v3;
}

// measure magnitude of provided vector
KOKKOS_FUNCTION void magnitude(CUDAREAL *vector) {
#ifdef __CUDA_ARCH__
        vector[0] = norm3d(vector[1], vector[2], vector[3]);
#else
        vector[0] = sqrt(vector[1] * vector[1] + vector[2] * vector[2] + vector[3] * vector[3]);
#endif
}

// scale magnitude of provided vector
KOKKOS_FUNCTION CUDAREAL vector_scale(CUDAREAL *vector, CUDAREAL *new_vector, CUDAREAL scale) {

        new_vector[1] = scale * vector[1];
        new_vector[2] = scale * vector[2];
        new_vector[3] = scale * vector[3];
        magnitude(new_vector);

        return new_vector[0];
}

// Fourier transform of a grating
KOKKOS_FUNCTION CUDAREAL sincg(CUDAREAL x, CUDAREAL N) {
        if (x != 0.0) {
                return sin(x * N) / sin(x);
        }
        return N;

}

#ifndef sinpi
#define  sinpi(x) sin(x * M_PI)
#endif

KOKKOS_INLINE_FUNCTION CUDAREAL sincgrad(CUDAREAL x, CUDAREAL N) {
        if (x != 0.0) {
                return sinpi(x * N) / sinpi(x);
        }
        return N;
}

// Fourier transform of a sphere
KOKKOS_FUNCTION CUDAREAL sinc3(CUDAREAL x) {
        if (x != 0.0) {
                return 3.0 * (sin(x) / x - cos(x)) / (x * x);
        }
        return 1.0;

}

#endif // SIMTBX_KOKKOS_KERNEL_MATH_H
