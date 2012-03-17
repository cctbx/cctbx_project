# A class for producing efficient looping limits for reflection
# prediction based on the Reeke algorithm (see Mosflm).

from scitbx import matrix
import scitbx.math
import math

class reeke_model:
    """Model and methods for the Reeke algorithm"""

    def __init__(self, ub, axis, s0, dmin, phi_beg, phi_end, margin = 3):

        # the original orientation, at phi = 0
        self._ub = ub

        # mapping of permuted axes p, q, and r
        self._permutation = None

        # the source vector and wavelength
        self._s0 = s0
        self._wavelength = 1 / math.sqrt(s0.dot(s0))

        # the rotation axis and angular range
        self._axis = axis
        self._phi_range = (phi_beg, phi_end)

        # the resolution limit
        self._dstarmax = 1 / dmin
        self._dstarmax2 = self._dstarmax**2

        # Margin by which to expand limits. Mosflm uses 3.
        # It might be useful to account for errors in the orientation.
        self._margin = int(margin)

        # Set the orientation at the beginning and end of this wedge, and
        # the rotation matrix for the beginning to the mid point.
        # NB The wedge could be the oscillation range for a single image, but
        # would need expanding by the rocking width and beam divergence to
        # catch all the partials on this image. Alternatively, expand only
        # at the extrema of the whole sweep and allow the partials to be
        # captured on adjacent wedges.

        r_beg = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
            axis = self._axis, angle = phi_beg, deg = True))
        self._r_half_osc = matrix.sqr(
            scitbx.math.r3_rotation_axis_and_angle_as_matrix(
            axis = self._axis, angle = (phi_end - phi_beg) / 2.0, deg=True))

        ub_beg = r_beg * self._ub
        ub_mid = self._r_half_osc * ub_beg
        ub_end = self._r_half_osc * ub_mid

        # Determine the permutation order of columns of the orientation
        # matrix. Use the orientation from the middle of the wedge for this.
        # As a side-effect set self._permutation.

        self._permutation = None
        col1, col2, col3 = self._permute_axes(ub_mid)

        # Thus set the reciprocal lattice axis vectors, in permuted order
        # p, q and r for both orientations

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

        # Set permuted orientation matrices

        self._p_beg = matrix.sqr(self._rlv_beg[0].elems +
                                 self._rlv_beg[1].elems +
                                 self._rlv_beg[2].elems).transpose()
        self._p_end = matrix.sqr(self._rlv_end[0].elems +
                                 self._rlv_end[1].elems +
                                 self._rlv_end[2].elems).transpose()

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

        # The following are set during the generation of indices

        # planes of constant p tangential to the Ewald sphere
        self._ewald_p_lim_beg = None
        self._ewald_p_lim_end = None

        # planes of constant p touching the circle of intersection between
        # the Ewald and resolution limiting spheres
        self._res_p_lim_beg = None
        self._res_p_lim_end = None

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

        # permutation matrix such that h, k, l = M * (p, q, r)
        elems = [int(0)] * 9
        elems[col1] = int(1)
        elems[col2 + 3] = int(1)
        elems[col3 + 6] = int(1)

        self._permutation = matrix.sqr(elems)

        # Return the permuted order of the columns

        return col1, col2, col3

    def _solve_quad(self, a, b, c):
        """Robust solution, for real roots only, of a quadratic in the form
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

        # TODO better way to get sin_2theta?
        sin_theta = 0.5 * self._wavelength * self._dstarmax
        assert abs(sin_theta) <= 1.0 # sanity check
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
        c =  cq[1] + q * (cq[2] + self._cp[13][0]) + \
             q**2 * self._cp[10] + cq[3][0]

        ewald_r_lim_beg = self._solve_quad(a, b, c)
        ewald_r_lim_beg = [item for item in ewald_r_lim_beg \
                           if item is not None]

        # Ewald sphere limits for the end orientation

        b = cq[0] + q * self._cp[8] + self._cp[12][0]
        c =  cq[1] + q * (cq[2] + self._cp[13][1]) + \
             q**2 * self._cp[10] + cq[3][1]

        ewald_r_lim_end = self._solve_quad(a, b, c)
        ewald_r_lim_end = [item for item in ewald_r_lim_end \
                           if item is not None]

        # if no intersections at all, return None

        if len(ewald_r_lim_beg) == 0 and len(ewald_r_lim_end) == 0:
            return None

        # if there are no intersections at the beginning orientation, set
        # a single loop covering the range between the intersections at
        # the end orientation, and vice versa.

        if len(ewald_r_lim_beg) == 0:

            l1 = [int(min(ewald_r_lim_end)) - max(self._margin, 1), \
                  int(max(ewald_r_lim_end)) + max(self._margin, 1)]
            l2 = [None]

        elif len(ewald_r_lim_end) == 0:

            l1 = [int(min(ewald_r_lim_beg)) - max(self._margin, 1), \
                  int(max(ewald_r_lim_beg)) + max(self._margin, 1)]
            l2 = [None]

        # otherwise there is at least one intersection at both orientations.
        # Set two loops, one for each range swept out by a point of
        # intersection as it travels from the beginning to the end
        # orientation.

        else:

            l1 = sorted([min(ewald_r_lim_beg), min(ewald_r_lim_end)])
            l1 = [int(l1[0]) - max(self._margin, 1), \
                  int(l1[1]) + max(self._margin, 1)]
            l2 = sorted([max(ewald_r_lim_beg), max(ewald_r_lim_end)])
            l2 = [int(l2[0]) - max(self._margin, 1), \
                  int(l2[1]) + max(self._margin, 1)]

        # restrict loops according to the resolution limit

        l1[0] = max(res_r_lim[0], l1[0])
        l1[1] = min(res_r_lim[1], l1[1])
        if l1[0] >= l1[1]: l1 = [None]

        if l2 != [None]:
            l2[0] = max(res_r_lim[0], l2[0])
            l2[1] = min(res_r_lim[1], l2[1])
            if l2[0] >= l2[1]: l2 = [None]

        if l1 == [None] and l2 == [None]: return None

        return [tuple(l1), tuple(l2)]

    def generate_indices(self):
        """Determine looping limits for indices h, k and l using the
        Reeke algorithm. This is the top level method for this module.
        All other methods are (probably) called by this, and therefore
        may as well be private. Can clean this up later. Also lots of
        debugging print statements to remove!"""

        # The outer loop is between limits for the axis most closely parallel,
        # or antiparallel, to the X-ray beam, which is called 'p'.

        # Determine the limiting values of p

        p_lim = self._p_limits()

        # fill indices list by looping over p, q and r

        hkl = []

        for p in range(p_lim[0], p_lim[1] + 1):

            # quantities that vary with p but are constant with q

            cq = [(p * self._cp[7]), \
                  (p**2 * self._cp[9]), \
                  (p * self._cp[11]), \
                  (p * self._cp[14][0], p * self._cp[14][1])]

            # find the limiting values of q

            q_lim = self._q_limits(p)
            if q_lim is None: continue

            for q in range(q_lim[0], q_lim[1] + 1):

                # find the limiting values of r

                r_lim = self._r_limits(p, q, cq)
                if r_lim is None: continue

                for item in r_lim:

                    if item[0] is None: continue

                    for r in range(item[0], item[1]+1):
                        hkl.append((self._permutation * (p, q, r)).elems)

        return hkl

    def visualize_with_rgl(self, rscript="reeke_vis.R", dat="reeke_hkl.dat"):
        """Write an R script and an associated data file
        for visualisation of generated indices between phi_beg and phi_end,
        using R and the rgl add-on package."""

        # Sorry, this is ugly. I don't know matplotlib yet.

        # write R script

        f = open(rscript, "w")

        f.write("# Run this from within R using\n" + \
                "# install.packages('rgl')\n" + \
                "# source('%s')\n\n" % rscript)
        f.write("library(rgl)\n")
        f.write("p_ax <- c(%.9f, %.9f, %.9f)\n" % self._rlv_beg[0].elems)
        f.write("q_ax <- c(%.9f, %.9f, %.9f)\n" % self._rlv_beg[1].elems)
        f.write("r_ax <- c(%.9f, %.9f, %.9f)\n" % self._rlv_beg[2].elems)
        f.write("s0 <- c(%.9f, %.9f, %.9f)\n" % self._s0.elems)
        f.write("rot_ax <- c(%.9f, %.9f, %.9f)\n" % self._axis.elems)
        f.write("phi_range <- c(%.9f, %.9f)\n" % self._phi_range)
        f.write("half_osc <- matrix(data = c(" + \
                "%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f,%.9f" % \
                self._r_half_osc.elems + \
                "), nrow=3, byrow=T)\n")
        f.write("dstarmax <- %.9f\n" % self._dstarmax)
        f.write("perm <- solve(matrix(data = c(" + \
                "%d,%d,%d,%d,%d,%d,%d,%d,%d" % self._permutation.elems + \
                "), nrow=3, byrow=T))\n")
        f.write("\n# draw the Ewald and limiting spheres\n" + \
                "open3d()\n" + \
                "spheres3d(s0,radius=sqrt(sum(s0*s0)),color='#CCCCFF'," + \
                "alpha=0.3)\n" + \
                "spheres3d(c(0,0,0),radius=dstarmax," + \
                "color='red',alpha=0.1)\n" + \
                "\n# draw the source vector and rotation axis\n" + \
                "lines3d(rbind(c(0,0,0),s0))\n" + \
                "lines3d(rbind(c(0,0,0),rot_ax))\n" + \
                "\n# draw the reciprocal lattice axes at ten times their " + \
                "length\n" + \
                "lines3d(rbind(c(0,0,0),10*p_ax),col='red')\n" + \
                "lines3d(rbind(c(0,0,0),10*q_ax),col='green')\n" + \
                "lines3d(rbind(c(0,0,0),10*r_ax),col='blue')\n"
                "s0mag <- sqrt(sum(s0*s0))\n" + \
                "s0unit <- s0 / s0mag\n" + \
                "\n# two unit vectors orthogonal to s0\n" + \
                "s0unitx1 <- c((-s0unit[3]),(0),(s0unit[1]))\n" + \
                "s0unitx2 <- c((s0unit[2]*s0unitx1[3] - " + \
                "s0unit[3]*s0unitx1[2]),\n" + \
                "   (s0unit[1]*s0unitx1[3] - s0unit[3]*s0unitx1[1]),\n" + \
                "   (s0unit[1]*s0unitx1[2] - s0unit[2]*s0unitx1[1]))\n" + \
                "sin_theta <- dstarmax/(2*s0mag)\n" + \
                "sin_2theta <- sin(2*asin(sin_theta))\n" + \
                "\n# distance to the centre of the circle of" + \
                "intersection, along s0\n" + \
                "e <- 2 * sqrt(sum(s0*s0)) * sin_theta ^2\n" + \
                "\n# radius of the circle of intersection\n" + \
                "R <- s0mag * sin_2theta\n" + \
                "\n# make points around the circle\n" + \
                "tau <- seq(from=0,to=2*pi,by=0.01)\n" + \
                "circ <- t(sapply(tau, function(x){\n" + \
                "  e * s0unit + R*sin(x) * s0unitx1 + R*cos(x) * s0unitx2}))\n" + \
                "\n# draw the circle\n" + \
                "lines3d(circ)\n" + \
                "\n# load the generated indices\n"
                "pts <- read.csv('./%s')\n" % dat)
        f.write("\n# convert h, k, l to reciprocal space coordinates\n" + \
                "conv <- function(h) {p <- perm %*% h\n" + \
                "    return(p[1]*p_ax + p[2]*q_ax + p[3]*r_ax)}\n" + \
                "pts <- t(apply(pts, MARGIN = 1, FUN = conv))\n" + \
                "\n# draw the generated indices\n" + \
                "points3d(pts, col='blue')\n" + \
                "\n")
        f.close()

        # write data file

        indices = self.generate_indices()

        f = open(dat, "w")
        f.write("h, k, l\n")

        for hkl in indices:
            f.write("%d, %d, %d\n" % hkl)
        f.close()

        print "Generated indices were written to %s" % dat
        print "An R script for visualising these was written to %s," % rscript
        print "which can be run from the R prompt with:"
        print "source('%s')" % rscript

        return

def reeke_model_for_use_case(phi_beg, phi_end, margin):
    """Construct a reeke_model for the geometry of the Use Case Thaumatin
    dataset, taken from the XDS XPARM. The values are hard-
    coded here so that this module does not rely on the location of that
    file."""

    axis = matrix.col([0.0, 1.0, 0.0])
    ub = matrix.sqr([-0.0133393674072, -0.00541609051856, -0.00367748834997,
                    0.00989309470346, 0.000574825936669, -0.0054505379664,
                    0.00475395109417, -0.0163935257377, 0.00102384915696])
    s0 = matrix.col([0.00237878589035, 1.55544539299e-16, -1.09015329696])
    dmin = 1.20117776325

    return reeke_model(ub, axis, s0, dmin, phi_beg, phi_end, margin)

def regression_test():
    """Perform a regression test by generating indices for a small wedge in
    phi for the Use Case data, and compare to expected results."""

    r = reeke_model_for_use_case(0, 0.2, 3)
    indices = r.generate_indices()

    h, k, l = zip(*indices)
    assert (min(h), max(h)) == (-47, 37)
    assert (min(k), max(k)) == (0, 32)
    assert (min(l), max(l)) == (-124, 109)

    # test reduced margin size
    r = reeke_model_for_use_case(0, 0.2, 1)
    indices = r.generate_indices()

    h, k, l = zip(*indices)
    assert (min(h), max(h)) == (-47, 37)
    assert (min(k), max(k)) == (2, 31)
    assert (min(l), max(l)) == (-122, 107)

    # test increased wedge angle
    r = reeke_model_for_use_case(0, 1.0, 3)
    indices = r.generate_indices()

    h, k, l = zip(*indices)
    assert (min(h), max(h)) == (-47, 37)
    assert (min(k), max(k)) == (-1, 32)
    assert (min(l), max(l)) == (-124, 109)

    #TODO Tests for an oblique cell
    #TODO Better tests than ranges of generated indices

    print "OK"

if __name__ == '__main__':

    import sys

    if len(sys.argv) == 1:
        regression_test()

    elif len(sys.argv) < 3:
        from libtbx.utils import Sorry
        raise Sorry("Expecting either 2 or 3 arguments: start_phi end_phi margin=3")

    else:

        phi_beg, phi_end = float(sys.argv[1]), float(sys.argv[2])
        margin = int(sys.argv[3]) if len(sys.argv) == 4 else 3

        # test run for development/debugging. Values come from the use case data
        r = reeke_model_for_use_case(phi_beg, phi_end, margin)
        indices = r.generate_indices()

        for hkl in indices:
            print "%4d %4d %4d" % hkl
