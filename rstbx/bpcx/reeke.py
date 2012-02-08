# A class for producing efficient looping limits for reflection
# prediction based on the Reeke algorithm (see Mosflm).

from scitbx import matrix
import scitbx.math
import math

class reeke_model:
    """Model and methods for the Reeke algorithm"""

    def __init__(self, ub, axis, s0):

        # the original orientation, at phi = 0
        self._ub = ub

        # mapping of permuted axes P, Q, and R
        self._permutation = 0, 1, 2

        # the source vector and wavelength
        self._s0 = s0

        # the rotation axis
        self._axis = axis

        # in here initialize other variables which you will subsequently
        # define and use - set them to None perhaps?

        # reciprocal lattice vectors at the beginning and end of the frame

        self._rlv_beg = None
        self._rlv_end = None

        return

    def get_s0(self):
        return self._s0

    def get_ub(self):
        return self._ub

    def get_axis(self):
        return self._axis

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

    def _ewald_p_limit(self):
        """Calculate the value of p at which planes of constant p are
        tangential to the Ewald sphere. Note p is the reciprocal cell
        axis given by the first column of the permuted orientation matrix."""

        v_beg = self._rlv_beg[1].cross(self._rlv_beg[2]).normalize()
        v_end = self._rlv_end[1].cross(self._rlv_end[2]).normalize()

        # Find distance between the planes of p

        p_dist = abs(self._rlv_beg[0].dot(v_beg))

        # Find distances between p = 0 and the centre of the Ewald sphere

        dp_beg = abs(v_beg.dot(self._s0))
        dp_end = abs(v_end.dot(self._s0))

        # There are two planes of constant p that are tangential to the Ewald
        # sphere, on either side of the sphere. The smaller in magnitude of p
        # is the number of planes that fit in one radius of the Ewald sphere
        # minus the number of planes between the centre of the Ewald sphere
        # and the p=0 plane (a diagram helps!). The larger is the number of
        # planes in one radius of the Ewald sphere *plus* the the number of
        # planes between the centre of the Ewald sphere and p=0.
        #
        # The correct sign is determined by whether p is more closely parallel
        # or antiparallel to the beam direction. If p increases along the beam
        # (i.e. against s0), then the smaller limit is positive.

        # The corner cases which we have discussed are mostly impossible
        # as we know rlv_beg[0] is close to colinear with the beam. This is
        # always going to be nasty.

        sign = cmp(self._rlv_beg[0].dot(self._s0), 0)

        limits = [(sign * s * (self._s0.length() + s *  dp_beg) / p_dist) \
                  for s in (-1, 1)]

        p_lim_beg = tuple(sorted(limits))

        sign = cmp(self._rlv_end[0].dot(self._s0), 0)

        limits = [(sign * s * (self._s0.length() + s *  dp_end) / p_dist) \
                  for s in (-1, 1)]

        p_lim_end = tuple(sorted(limits))

        return p_lim_beg, p_lim_end

    # This method will probably return the list of indices, not the limits, in
    # which case it should be renamed 'generate_indices_reeke', or similar.

    def loop_limits(self, phi_beg, phi_end):
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

        # perhaps a better convention is:
        #
        # print 'parameter "%s" is nonsense' % parameter

        print "The UB matrix after rotation to the centre of the " + \
              "wedge at %.3f degrees is" % \
              round(phi_beg + (phi_end - phi_beg)/2.0, 5)
        print ub_mid.round(5)
        print "The reciprocal cell axes are permuted in order %d %d %d" % \
              self._permutation
        print "giving a unit vector in the p direction,"

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
        self.p_beg = matrix.sqr(self._rlv_beg[0].elems +
                                self._rlv_beg[1].elems +
                                self._rlv_beg[2].elems).transpose()
        self.p_end = matrix.sqr(self._rlv_end[0].elems +
                                self._rlv_end[1].elems +
                                self._rlv_end[2].elems).transpose()
        #self.p = p_beg, p_end

        print "The permuted orientation matrix at the beginning of the wedge is"
        print self.p_beg.round(5)
        print "and the permuted orientation matrix at the end of the wedge is"
        print self.p_end.round(5)
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
        pp_beg = matrix.rec(self.p_beg.elems[0:3] + (-1.*self._s0[0],) +
                              self.p_beg.elems[3:6] + (-1.*self._s0[1],) +
                              self.p_beg.elems[6:9] + (-1.*self._s0[2],), n=(3, 4))
        pp_end = matrix.rec(self.p_end.elems[0:3] + (-1.*self._s0[0],) +
                            self.p_end.elems[3:6] + (-1.*self._s0[1],) +
                            self.p_end.elems[6:9] + (-1.*self._s0[2],), n=(3, 4))
        #self.pp = pp_beg, pp_end

        # Various quantities of interest are obtained from the reciprocal metric
        # tensor T of P'. These quantities are to be used (later) for solving the
        # intersection of a line of constant p, q index with the Ewald sphere. It
        # is efficient to calculate these before the outer loop. So, calculate T
        # for both beginning and end orientations
        #self.t = map(lambda m: m.transpose() * m, self.pp)
        t_beg = pp_beg.transpose() * pp_beg
        t_end = pp_end.transpose() * pp_end

        # The outer loop is between limits for the axis most closely parallel,
        # or antiparallel, to the X-ray beam, which is called 'p'.
        print "Ewald sphere limits for p index are p ="
        p_lim = self._ewald_p_limit()
        print "%.3f, %.3f for the start orientation" % p_lim[0]
        print "%.3f, %.3f for the end orientation" % p_lim[1]

        print "The axis in the p direction is"
        print self._rlv_beg[0].round(5)
        print "and the vector s0 is"
        print self._s0.round(5)
        print "so the unit vector in the beam direction is"
        print -1.0 * self._s0.normalize().round(5)

        return


if __name__ == '__main__':

    # test run for development/debugging.
    axis = matrix.col([0.0, 1.0, 0.0])
    ub = matrix.sqr([-0.0133393674072, -0.00541609051856, -0.00367748834997,
                    0.00989309470346, 0.000574825936669, -0.0054505379664,
                    0.00475395109417, -0.0163935257377, 0.00102384915696])
    s0 = matrix.col([0.00237878589035, 1.55544539299e-16,-1.09015329696])
    r = reeke_model(ub, axis, s0)

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

    r.loop_limits(30, 30.1)
