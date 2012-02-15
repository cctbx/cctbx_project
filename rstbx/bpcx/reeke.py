# A class for producing efficient looping limits for reflection
# prediction based on the Reeke algorithm (see Mosflm).

from scitbx import matrix
import scitbx.math
import math

class reeke_model:
    """Model and methods for the Reeke algorithm"""

    def __init__(self, ub, axis, s0, dmin, margin = 3):

        # the original orientation, at phi = 0
        self._ub = ub

        # mapping of permuted axes p, q, and r
        self._permutation = 0, 1, 2

        # the source vector and wavelength
        self._s0 = s0
        self._wavelength = 1 / math.sqrt(s0.dot(s0))

        # the rotation axis
        self._axis = axis

        # the resolution limit
        self._dstarmax = 1 / dmin
        self._dstarmax2 = self._dstarmax**2

        # margin by which to expand limits. Mosflm uses 3. Can play with
        # this value to see what difference it makes - perhaps it is to
        # account for errors in the orientation
        self._margin = int(margin)

        # reciprocal lattice vectors at the beginning and end of the wedge

        self._rlv_beg = None
        self._rlv_end = None

        # orientation at beginning and end of wedge, in permuted order
        self._p_beg = None
        self._p_end = None

        # planes of constant p tangential to the Ewald sphere

        self._ewald_p_lim_beg = None
        self._ewald_p_lim_end = None

        # planes of constant p touching the circle of intersection between
        # the Ewald and resolution limiting spheres

        self._res_p_lim_beg = None
        self._res_p_lim_end = None

        # quantities that are constant with p
        self._cp = None

        # looping p limits

        self._p_lim = None

        return

    def get_s0(self):
        return self._s0

    def get_ub(self):
        return self._ub

    def get_axis(self):
        return self._axis

    def get_all_p_limits(self):
        """Get both the Ewald and limiting sphere limits for planes of p.
        This is useful for plotting the planes, for example."""
        return (self._ewald_p_lim_beg, self._ewald_p_lim_end, \
                self._res_p_lim_beg, self._res_p_lim_end)

    def _permute_axes(self, ub):
        """Find permutation of the columns of an orientation matrix so that
        column p is closest to the X-ray beam direction s0, column r is
        closest of q and r to the spindle axis and column q is the remaining
        direction."""

        # Extract the reciprocal lattice directions from the columns of UB

        rl_dirs = [matrix.col(v).normalize() for v in \
                   ub.transpose().as_list_of_lists()]

        # Find reciprocal lattice axis closest to s0 by checking magnitude
        # of dot products between normalised axes and s0, then swap as required

        along_beam = [math.fabs(rl_dirs[j].dot(self._s0)) for j in range(3)]

        col1 = along_beam.index(max(along_beam))

        rl_dirs[0], rl_dirs[col1] = rl_dirs[col1], rl_dirs[0]

        # Now find which of the two remaining reciprocal lattice axes is
        # closest to the rotation axis.

        along_spindle = [math.fabs(rl_dirs[j].dot(self._axis)) for j in (1, 2)]

        col3 = along_spindle.index(max(along_spindle)) + 1

        # Which is the remaining column index?

        col2 = [j for j in range(3) if not j in (col1, col3)][0]

        # couldn't this be stored as a matrix which would mean you could
        # have h, k, l = M * (p, q, r)

        self._permutation = col1, col2, col3

        # Return the permuted order of the columns

        return col1, col2, col3

    def _solve_quad(self, a, b, c):
        """Robust solution for real roots only of a quadratic in the form
        (ax^2 + bx + c)."""

        discriminant = b**2 - 4 * a * c

        if discriminant > 0:
            sign = cmp(b, 0)
            if sign == 0: sign = 1.0
            q = -0.5 * (b + sign * math.sqrt(discriminant))
            x1 = q / a if a != 0 else None
            x2 = c / q if q != 0 else None
            return [x1, x2]

        elif discriminant == 0:
            return [(-b) / (2 * a)] * 2

        else:
            return [None]


    def _p_limits(self):
        """Calculate the values of p at which planes of constant p are
        tangential to the Ewald sphere, and values of p at which planes
        of constant p touch the circle of intersection between the Ewald
        and resolution limiting sphere. Note p is the reciprocal cell
        axis given by the first column of the permuted orientation matrix.
        Set the limits as attributes and return a single set of overall
        limits."""

        # Calculate unit vectors normal to planes of constant p, ensuring
        # they point in the direction of increasing p.

        v_beg = self._rlv_beg[1].cross(self._rlv_beg[2]).normalize()

        if self._rlv_beg[0].dot(v_beg) < 0:
            v_beg = -1 * v_beg

        v_end = self._rlv_end[1].cross(self._rlv_end[2]).normalize()

        if self._rlv_end[0].dot(v_end) < 0:
            v_end = -1 * v_end

        # Find distance between the planes of p

        p_dist = abs(self._rlv_beg[0].dot(v_beg))

        # Find distances between p = 0 and the plane passing through the
        # centre of the Ewald sphere

        dp_beg = abs(v_beg.dot(self._s0))
        dp_end = abs(v_end.dot(self._s0))

        # There are two planes of constant p that are tangential to the Ewald
        # sphere, on either side of the sphere. The smaller in magnitude of p
        # is the number of planes that fit in one radius of the Ewald sphere
        # minus the number of planes between the centre of the Ewald sphere
        # and the p=0 plane (a diagram helps!). The larger is the number of
        # planes in one radius of the Ewald sphere *plus* the the number of
        # planes between the centre of the Ewald sphere and p = 0.
        #
        # The correct sign is determined by whether the plane normal vector is
        # more closely parallel or antiparallel to the beam direction.

        sign = cmp(v_beg.dot(self._s0), 0)

        limits = [(sign * s * (self._s0.length() + s *  dp_beg) / p_dist) \
                  for s in (-1, 1)]

        self._ewald_p_lim_beg = tuple(sorted(limits))

        sign = cmp(v_end.dot(self._s0), 0)

        limits = [(sign * s * (self._s0.length() + s *  dp_end) / p_dist) \
                  for s in (-1, 1)]

        self._ewald_p_lim_end = tuple(sorted(limits))

        # Now determine limits for the planes of p that touch the circle of
        # intersection between the Ewald and resolution limiting spheres

        # better way to get sin_2theta?
        # also need some sanity checks in case dmin is ridiculous
        sin_theta = 0.5 * self._wavelength * self._dstarmax
        sin_2theta = math.sin(2.0 * math.asin(sin_theta))

        e = 2.0 * sin_theta**2 * dp_beg
        f = sin_2theta * math.sqrt(max(1.0 / self._wavelength**2 - dp_beg**2, 0))
        limits = [(sign * e + s * f) / p_dist for s in (-1, 1)]

        self._res_p_lim_beg = tuple(sorted(limits))

        e = 2.0 * sin_theta**2 * dp_end
        f = sin_2theta * math.sqrt(max(1.0 / self._wavelength**2 - dp_end**2, 0))
        limits = [(sign * e + s * f) / p_dist for s in (-1, 1)]

        self._res_p_lim_end = tuple(sorted(limits))

        # select most restrictive of Ewald and resolution limits
        p_lim_beg = sorted(self._ewald_p_lim_beg + self._res_p_lim_beg)[1:3]
        p_lim_end = sorted(self._ewald_p_lim_end + self._res_p_lim_end)[1:3]

        # single set of limits covering overall range
        p_lim = sorted(p_lim_beg + p_lim_end)[0::3]
        p_lim[0] = int(p_lim[0]) - self._margin
        p_lim[1] = int(p_lim[1]) + self._margin

        return p_lim

    def _q_limits(self, p):
        """Calculate the values of q at which lines of constant p, q are
        tangential to the circle intersecting the Ewald sphere at plane p,
        and values of q at which lines of constant p, q are tangential to
        the circle intersecting the resolution limiting sphere at plane p.i
        Return the appropriate overall limits."""

        # First the resolution limits. Set up the quadratic to solve

        a = self._cp[6]
        b = 2.0 * p * self._cp[5]
        c = p**2 * self._cp[3] + self._cp[0] * self._dstarmax2

        res_q_lim = self._solve_quad(a, b, c)
        res_q_lim = sorted([item for item in res_q_lim \
                            if item is not None])
        if len(res_q_lim) == 0: return None

        # Extend limits by the margin, ensuring there is a range even for
        # a single quadratic root
        res_q_lim = [int(res_q_lim[0]) - max(self._margin, 1),
                     int(res_q_lim[-1]) + max(self._margin, 1)]

        # Ewald sphere limits for the beginning orientation

        b = 2.0 * (self._cp[4][0] + p * self._cp[5])
        c = self._cp[1][0] + p * (2 * self._cp[2][0] + p * self._cp[3])

        ewald_q_lim_beg = self._solve_quad(a, b, c)

        # Ewald sphere limits for the end orientation

        b = 2.0 * (self._cp[4][1] + p * self._cp[5])
        c = self._cp[1][1] + p * (2 * self._cp[2][1] + p * self._cp[3])

        ewald_q_lim_end = self._solve_quad(a, b, c)

        # Determine the overall Ewald limits
        ewald_q_lim = sorted([item for item in ewald_q_lim_beg + \
                              ewald_q_lim_end if item is not None])
        if len(ewald_q_lim) > 0:
            ewald_q_lim = [int(ewald_q_lim[0]) - max(self._margin, 1),
                           int(ewald_q_lim[-1]) + max(self._margin, 1)]

        else:
            return None

        # Choose most restrictive of Ewald and res limits. The expansion of
        # limits by the margin ensures that we have a 4 element list here

        q_lim = sorted(res_q_lim + ewald_q_lim)
        q_lim = [q_lim[1], q_lim[2]]

        return q_lim

    def _r_limits(self, p, q, cq):
        """Calculate the values of r at which lines of constant p, q
        intersect the resolution limiting and the Ewald spheres, and
        return the appropriate overall limits"""

        # First the resolution limits. Set up the quadratic to solve

        a = self._cp[0]
        b = cq[0] + q * self._cp[8]
        c = cq[1] + q**2 * self._cp[10] + q * cq[2] - self._dstarmax2

        res_r_lim = self._solve_quad(a, b, c)
        res_r_lim = sorted([item for item in res_r_lim if item is not None])
        if len(res_r_lim) == 0: return None

        # Extend limits by the margin, ensuring there is a range even for
        # a single quadratic root
        res_r_lim = [int(res_r_lim[0]) - max(self._margin, 1),
                     int(res_r_lim[-1]) + max(self._margin, 1)]

       # Ewald sphere limits for the beginning orientation

        b =  cq[0] + q * self._cp[8] + self._cp[12][0]
        c =  cq[1] + q * (cq[2] + self._cp[13][0]) + q**2 * self._cp[10] + cq[3][0]

        ewald_r_lim_beg = self._solve_quad(a, b, c)

        # Alternative way to calculate - as used by Mosflm. Are they equal?

        # Ewald sphere limits for the end orientation

        b = cq[0] + q * self._cp[8] + self._cp[12][0]
        c =  cq[1] + q * (cq[2] + self._cp[13][1]) + q**2 * self._cp[10] + cq[3][1]

        ewald_r_lim_end = self._solve_quad(a, b, c)

        return ewald_r_lim_beg, ewald_r_lim_end

    def generate_indices(self, phi_beg, phi_end):
        """Determine looping limits for indices h, k and l using the
        Reeke algorithm. This is the top level method for this module.
        All other methods are (probably) called by this, and therefore
        may as well be private. Can clean this up later. Also lots of
        debugging print statements to remove!"""

        # First set the orientation at the beginning and end of this wedge
        # (typically this is the oscillation range of a single image)
        # NB Mosflm extends the rotation range at each end by the maximum
        # reflection width, to catch all partials for this image

        r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
            axis = self._axis, angle = phi_beg, deg = True))
        r_half_osc = matrix.sqr(
            scitbx.math.r3_rotation_axis_and_angle_as_matrix(
            axis = self._axis, angle = (phi_end - phi_beg) / 2.0, deg=True))

        ub_beg = r_beg * self._ub
        ub_mid = r_half_osc * ub_beg
        ub_end = r_half_osc * ub_mid

        # Determine the permutation order of columns of the orientation matrix.
        # Use the orientation from the middle of the wedge for this.

        col1, col2, col3 = self._permute_axes(ub_mid)

        # Thus set the reciprocal lattice axis vectors, in permuted order p, q and r
        rl_vec = [ub_beg.extract_block(start=(0,0), stop=(3,1)),
                  ub_beg.extract_block(start=(0,1), stop=(3,2)),
                  ub_beg.extract_block(start=(0,2), stop=(3,3))]
        self._rlv_beg = [rl_vec[col1],
                         rl_vec[col2],
                         rl_vec[col3]]
        rl_vec = [ub_end.extract_block(start=(0,0), stop=(3,1)),
                  ub_end.extract_block(start=(0,1), stop=(3,2)),
                  ub_end.extract_block(start=(0,2), stop=(3,3))]
        self._rlv_end = [rl_vec[col1],
                         rl_vec[col2],
                         rl_vec[col3]]

        # Set permuted orientation matrices, for beginning and end of wedge
        self._p_beg = matrix.sqr(self._rlv_beg[0].elems +
                                self._rlv_beg[1].elems +
                                self._rlv_beg[2].elems).transpose()
        self._p_end = matrix.sqr(self._rlv_end[0].elems +
                                self._rlv_end[1].elems +
                                self._rlv_end[2].elems).transpose()
        #self.p = _p_beg, _p_end

        print "The permuted orientation matrix at the beginning of the wedge is"
        print self._p_beg.round(5)
        print "and the permuted orientation matrix at the end of the wedge is"
        print self._p_end.round(5)
        print

        # Define a new coordinate system concentric with the Ewald sphere.
        #
        # X' = X - s0_x
        # Y' = Y - s0_y
        # Z' = Z - s0_z
        #
        # X = P' h'
        # -   =  -
        #                                    p11 p12 p13 -s0_X
        # where h' = (p, q, r, 1)^T and P' = p21 p22 p23 -s0_y
        #       -                       =    p31 p32 p33 -s0_z
        #
        # Calculate P' matrices for the beginning and end orientations
        pp_beg = matrix.rec(self._p_beg.elems[0:3] + (-1.*self._s0[0],) +
                            self._p_beg.elems[3:6] + (-1.*self._s0[1],) +
                            self._p_beg.elems[6:9] + (-1.*self._s0[2],), n=(3, 4))
        pp_end = matrix.rec(self._p_end.elems[0:3] + (-1.*self._s0[0],) +
                            self._p_end.elems[3:6] + (-1.*self._s0[1],) +
                            self._p_end.elems[6:9] + (-1.*self._s0[2],), n=(3, 4))
        #self.pp = pp_beg, pp_end

        # Various quantities of interest are obtained from the reciprocal metric
        # tensor T of P'. These quantities are to be used (later) for solving the
        # intersection of a line of constant p, q index with the Ewald sphere. It
        # is efficient to calculate these before the outer loop. So, calculate T
        # for both beginning and end orientations
        t_beg = (pp_beg.transpose() * pp_beg).as_list_of_lists()
        t_end = (pp_end.transpose() * pp_end).as_list_of_lists()

        # quantities that are constant with p
        self._cp = [(t_beg[2][2]), \
                    (t_beg[2][3]**2, t_end[2][3]**2), \
                    (t_beg[0][2] * t_beg[2][3] - t_beg[0][3] * t_beg[2][2], \
                     t_end[0][2] * t_end[2][3] - t_end[0][3] * t_end[2][2]), \
                    (t_beg[0][2]**2 - t_beg[0][0] * t_beg[2][2]), \
                    (t_beg[1][2] * t_beg[2][3] - t_beg[1][3] * t_beg[2][2], \
                     t_end[1][2] * t_end[2][3] - t_end[1][3] * t_end[2][2]), \
                    (t_beg[0][2] * t_beg[1][2] - t_beg[0][1] * t_beg[2][2]),
                    (t_beg[1][2]**2 - t_beg[1][1] * t_beg[2][2]), \
                    (2.0 * t_beg[0][2]), \
                    (2.0 * t_beg[1][2]), \
                    (t_beg[0][0]), \
                    (t_beg[1][1]), \
                    (2.0 * t_beg[0][1]), \
                    (2.0 * t_beg[2][3], 2.0 * t_end[2][3]), \
                    (2.0 * t_beg[1][3], 2.0 * t_end[1][3]), \
                    (2.0 * t_beg[0][3], 2.0 * t_end[0][3])]
        #_cp[0]   t33
        #_cp[1]   R00_beg, R00_end
        #_cp[2]   R01_beg, R01_end
        #_cp[3]   R02
        #_cp[4]   R10_beg, R10_end
        #_cp[5]   R11
        #_cp[6]   R2
        #_cp[7]   2 t13
        #_cp[8]   2 t23
        #_cp[9]   t11
        #_cp[10]  t22
        #_cp[11]  2 t12
        #_cp[12]  2 t34_beg, 2 t34_end
        #_cp[13]  2 t24_beg, 2 t24_end
        #_cp[14]  2 t14_beg, 2 t14_end

        # The outer loop is between limits for the axis most closely parallel,
        # or antiparallel, to the X-ray beam, which is called 'p'.

        # Determine the limiting values of p
        p_lim = self._p_limits()
        x = self.get_all_p_limits()
        print "Ewald sphere limits for p index are p ="
        print "%.3f, %.3f for the start orientation" % x[0]
        print "%.3f, %.3f for the end orientation" % x[1]

        print "Limiting sphere limits for p index are p ="
        print "%.3f, %.3f for the start orientation" % x[2]
        print "%.3f, %.3f for the end orientation" % x[3]

        print "chosen looping limits", p_lim

        # loop over p
        for p in range(p_lim[0], p_lim[1] + 1):

            # quantities that vary with p but are constant with q
            cq = [(p * self._cp[7]), \
                  (p**2 * self._cp[9]), \
                  (p * self._cp[11]), \
                  (p * self._cp[14][0], p * self._cp[14][1])]
            #cq[0]    2 p t13
            #cq[1]    p^2 t11
            #cq[2]    2 p t12
            #cq[3]    2 p t14

            # find the limiting values of q
            q_lim = self._q_limits(p)

            # loop over q
            if q_lim is None: continue
            for q in range(q_lim[0], q_lim[1] + 1):

                # find the limiting values of r
                r_lim = self._r_limits(p, q, cq)
                print "for p = %d, q = %d, r_lim are " % (p, q), r_lim
        return


if __name__ == '__main__':

    # test run for development/debugging.
    axis = matrix.col([0.0, 1.0, 0.0])
    ub = matrix.sqr([-0.0133393674072, -0.00541609051856, -0.00367748834997,
                    0.00989309470346, 0.000574825936669, -0.0054505379664,
                    0.00475395109417, -0.0163935257377, 0.00102384915696])
    s0 = matrix.col([0.00237878589035, 1.55544539299e-16, 1.09015329696])
    dmin = 1.20117776325
    r = reeke_model(ub, axis, s0, dmin)

    print
    print "SETTING UP THE REEKE MODEL"
    print "s0 is"
    print r._s0.round(5)
    print "Original phi=0 orientation matrix UB is"
    print r._ub.round(5)
    print "Rotation axis is"
    print r._axis.round(5)
    print
    print "DETERMINING LOOP LIMITS FOR ROTATION BETWEEN 30 and 30.1 DEGREES"

    r.generate_indices(30, 30.1)
    #r.generate_indices(35, 35.1)
