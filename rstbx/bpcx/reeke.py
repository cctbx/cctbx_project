# A class for producing efficient looping limits for reflection
# prediction based on the Reeke algorithm (see Mosflm).

from scitbx import matrix
import scitbx.math

class reeke_model:
    """Model and methods for the Reeke algorithm"""

    def __init__(self, ub, axis, s0):

        # the original orientation, at phi = 0
        self.ub = ub

        # mapping of permuted axes P, Q, and R
        self.permutation = dict(P="h", Q="k", R="l")

        # the source vector and wavelength
        self.s0 = s0
        # wavelength is not required - take it from s0
        #self.wavelength = wavelength

        # the rotation axis
        self.axis = axis

    def _permute_axes(self, ub):
        """
        Find permutation of the columns of an orientation matrix so that
        column p is closest to the X-ray beam direction s0, column r is
        closest of q and r to the spindle axis and column q is the remaining
        direction
        """

        # Extract the reciprocal lattice direction vectors from the columns of UB
        rl_directions = [ub.extract_block(start=(0,0), stop=(3,1)).normalize(),
                         ub.extract_block(start=(0,1), stop=(3,2)).normalize(),
                         ub.extract_block(start=(0,2), stop=(3,3)).normalize()]

        # Find reciprocal lattice axis closest to s0 by checking magnitude
        # of dot products between normalised axes and s0, then reorder the list
        # as required
        along_beam = (abs(rl_directions[0].dot(self.s0)),
                      abs(rl_directions[1].dot(self.s0)),
                      abs(rl_directions[2].dot(self.s0)))
        col1 = along_beam.index(max(along_beam))
        rl_directions[0], rl_directions[col1] = rl_directions[col1], rl_directions[0]

        # Now find which of the two remaining reciprocal lattice axes is closest
        # to the rotation axis.
        along_spindle = (abs(rl_directions[1].dot(self.axis)),
                         abs(rl_directions[2].dot(self.axis)))
        col3 = along_spindle.index(max(along_spindle)) + 1

        # Keeping the reciprocal lattice unit directions might be useful? If not, can
        # remove the next two lines once the code is more complete.
        rl_directions[2], rl_directions[col3] = rl_directions[col3], rl_directions[2]
        self.rl_directions = rl_directions

        # Which is the remaining column index?
        col2 = [0, 1, 2]
        col2[col1] = 0
        col2[col3] = 0
        col2 = sum(col2)

        # Return the permuted order of the columns
        perm = ("h", "k", "l")
        self.permutation["P"] = perm[col1]
        self.permutation["Q"] = perm[col2]
        self.permutation["R"] = perm[col3]
        return col1, col2, col3

    def _ewald_p_limit(self):
        """
        Calculate the value of p at which planes of constant p are
        tangential to the Ewald sphere. Note p is the reciprocal cell
        axis given by the first column of the permuted orientation matrix
        """

        # Determine vector normal to the plane p = 0, by the cross product
        # of vectors q and r.
        #v = map(lambda v: v[1].cross(v[2]), self.rl_vec)
        v_start = self.rl_vectors_start[1].cross(self.rl_vectors_start[2]).normalize()
        v_end = self.rl_vectors_end[1].cross(self.rl_vectors_end[2]).normalize()

        # Find distance between the planes of p (this is the same for either orientation
        # so just do it for the start)
        p_dist = abs(self.rl_vectors_start[0].dot(v_start))

        # Find distance between plane of p = 0 and the centre of the Ewald sphere
        dp_start = abs(v_start.dot(self.s0))
        dp_end = abs(v_end.dot(self.s0))

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

        p_sign_start = 1.0 if self.rl_vectors_start[0].dot(self.s0) < 0 else -1.0
        p_lim1 = p_sign_start * (self.s0.length() - dp_start) / p_dist
        p_lim2 = -1.0 * p_sign_start * (self.s0.length() + dp_start) / p_dist
        p_lim_start = tuple(sorted([p_lim1, p_lim2]))

        p_sign_end = 1.0 if self.rl_vectors_end[0].dot(self.s0) < 0 else -1.0
        p_lim1 = p_sign_end * (self.s0.length() - dp_end) / p_dist
        p_lim2 = -1.0 * p_sign_end * (self.s0.length() + dp_end) / p_dist
        p_lim_end = tuple(sorted([p_lim1, p_lim2]))

        return p_lim_start, p_lim_end

    # This method will probably return the list of indices, not the limits, in
    # which case it should be renamed 'generate_indices_reeke', or similar.
    def loop_limits(self, phi_start, phi_end):
        """
        Determine looping limits for indices h, k and l using the
        Reeke algorithm. This is the top level method for this module.
        All other methods are (probably) called by this, and therefore
        may as well be private. Can clean this up later. Also lots of
        debugging print statements to remove!
        """

        # First set the orientation at the start and end of this wedge
        # (typically this is the oscillation range of a single image)
        # NB Mosflm extends the rotation range at each end by the maximum
        # reflection width, to catch all partials for this image

        r_start = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                axis=self.axis, angle=phi_start, deg=True))
        r_half_osc = matrix.sqr(scitbx.math.r3_rotation_axis_and_angle_as_matrix(
                axis=self.axis, angle=(phi_end - phi_start)/2.0, deg=True))

        ub_start = r_start * self.ub
        ub_mid = r_half_osc * ub_start
        ub_end = r_half_osc * ub_mid

        # Determine the permutation order of columns of the orientation matrix.
        # Use the orientation from the middle of the wedge for this.
        col1, col2, col3 = self._permute_axes(ub_mid)

        print "The UB matrix after rotation to the centre of the wedge at %.3f degrees is" % \
               round(phi_start + (phi_end - phi_start)/2.0, 5)
        print ub_mid.round(5)
        print ("The reciprocal cell axes are permuted in order " +
               self.permutation["P"] +
               self.permutation["Q"] +
               self.permutation["R"])
        print "giving a unit vector in the p direction,"
        print self.rl_directions[0].round(5)
        print "that is most nearly aligned with the beam direction,"
        print "a unit vector in the r direction,"
        print self.rl_directions[2].round(5)
        print "that is most nearly aligned with the spindle,"
        print "and a remaining unit vector in the q direction"
        print self.rl_directions[1].round(5)
        print

        # Thus set the reciprocal lattice axis vectors, in permuted order p, q and r
        rl_vec = [ub_start.extract_block(start=(0,0), stop=(3,1)),
                  ub_start.extract_block(start=(0,1), stop=(3,2)),
                  ub_start.extract_block(start=(0,2), stop=(3,3))]
        self.rl_vectors_start = [rl_vec[col1],
                                 rl_vec[col2],
                                 rl_vec[col3]]
        rl_vec = [ub_end.extract_block(start=(0,0), stop=(3,1)),
                  ub_end.extract_block(start=(0,1), stop=(3,2)),
                  ub_end.extract_block(start=(0,2), stop=(3,3))]
        self.rl_vectors_end = [rl_vec[col1],
                               rl_vec[col2],
                               rl_vec[col3]]

        # Set permuted orientation matrices, for beginning and end of wedge
        self.p_start = matrix.sqr(self.rl_vectors_start[0].elems +
                                  self.rl_vectors_start[1].elems +
                                  self.rl_vectors_start[2].elems).transpose()
        self.p_end = matrix.sqr(self.rl_vectors_end[0].elems +
                                self.rl_vectors_end[1].elems +
                                self.rl_vectors_end[2].elems).transpose()
        #self.p = p_start, p_end

        print "The permuted orientation matrix at the beginning of the wedge is"
        print self.p_start.round(5)
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
        # Calculate P' matrices for the start and end orientations
        pp_start = matrix.rec(self.p_start.elems[0:3] + (-1.*self.s0[0],) +
                              self.p_start.elems[3:6] + (-1.*self.s0[1],) +
                              self.p_start.elems[6:9] + (-1.*self.s0[2],), n=(3, 4))
        pp_end = matrix.rec(self.p_end.elems[0:3] + (-1.*self.s0[0],) +
                            self.p_end.elems[3:6] + (-1.*self.s0[1],) +
                            self.p_end.elems[6:9] + (-1.*self.s0[2],), n=(3, 4))
        #self.pp = pp_start, pp_end

        # Various quantities of interest are obtained from the reciprocal metric
        # tensor T of P'. These quantities are to be used (later) for solving the
        # intersection of a line of constant p, q index with the Ewald sphere. It
        # is efficient to calculate these before the outer loop. So, calculate T
        # for both start and end orientations
        #self.t = map(lambda m: m.transpose() * m, self.pp)
        t_start = pp_start.transpose() * pp_start
        t_end = pp_end.transpose() * pp_end

        # The outer loop is between limits for the axis most closely parallel,
        # or antiparallel, to the X-ray beam, which is called 'p'.
        print "Ewald sphere limits for p index are p ="
        p_lim = self._ewald_p_limit()
        print "%.3f, %.3f for the start orientation" % p_lim[0]
        print "%.3f, %.3f for the end orientation" % p_lim[1]

        print "The axis in the p direction is"
        print self.rl_vectors_start[0].round(5)
        print "and the vector s0 is"
        print self.s0.round(5)
        print "so the unit vector in the beam direction is"
        print -1.0 * self.s0.normalize().round(5)

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
    print r.s0.round(5)
    print "Original phi=0 orientation matrix UB is"
    print r.ub.round(5)
    print "Rotation axis is"
    print r.axis.round(5)
    print
    print "DETERMINING LOOP LIMITS FOR ROTATION BETWEEN 30 and 30.1 DEGREES"

    r.loop_limits(30, 30.1)
