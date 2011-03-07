#!/usr/bin/env python
import math

from scitbx import matrix
from rstbx.diffraction import partial_spot_position_partial_H
from libtbx.test_utils import approx_equal
from cctbx.crystal_orientation import crystal_orientation


some_indices = ([-5,29,0],[-5,20,0],[-3,13,0],[7,20,5],[10,10,10],[2,13,-5])

finite_differences = ([0.1,0.,0.],[0.,0.1,0.],[0.,0.,0.1])

# Reproduce the test provided by Graeme Winter)
def rotation_scattering(ra,rotation_vector,Amat,wavelength):
    beam_vector = matrix.col([0, 0, 1 / wavelength])

    for hkl in some_indices:
        if ra(hkl):

            omegas = ra.get_intersection_angles()

            for omegaidx in [0,1]:

                rot_mat = rotation_vector.axis_and_angle_as_r3_rotation_matrix(omegas[omegaidx])
                assert(-0.0001 < math.fabs(rot_mat.determinant() - 1.0) < 0.0001)

                H1 = (rot_mat * matrix.sqr(Amat))*matrix.col(hkl)
                H1 =  H1 + beam_vector
                len_H1 = math.sqrt((H1[0] * H1[0]) +
                                   (H1[1] * H1[1]) +
                                   (H1[2] * H1[2]))

                #print 'H1 wav %.3f %.3f' % (len_H1, 1.0 / wavelength)

                if math.fabs(len_H1 - 1.0 / wavelength) > 0.0001:
                    raise RuntimeError, 'length error for %d %d %d' % hkl

def test_finite_differences(ra,rotation_vector,Amat,wavelength):
    beam_vector = matrix.col([0, 0, 1 / wavelength])

    for hkl in [matrix.col(i) for i in some_indices]:
      if ra(hkl):
        omegas = ra.get_intersection_angles()
        dangles = [ra.dangle_(0),ra.dangle_(1)]

        for n,difference in enumerate([matrix.col(d) for d in finite_differences]):
          assert ra(hkl+difference) # not always true, true in these cases
          domegas = ra.get_intersection_angles()

          for omegaidx in [0,1]:
            #print "HKL","%18s"%(str(hkl.elems)),"idx",omegaidx, "%7.2f"%(180.*omegas[omegaidx]/math.pi),

            rot_mat = rotation_vector.axis_and_angle_as_r3_rotation_matrix(omegas[omegaidx])
            H1 = (rot_mat * matrix.sqr(Amat))*matrix.col(hkl) + beam_vector

            #print "%10.7f %10.7f %10.7f"%H1.elems
            #print "HKL","%18s"%str((hkl+difference).elems),"idx",omegaidx, "%7.2f"%(180.*domegas[omegaidx]/math.pi),

            rot_mat = rotation_vector.axis_and_angle_as_r3_rotation_matrix(domegas[omegaidx])
            H1 = (rot_mat * matrix.sqr(Amat))*matrix.col(hkl+difference) + beam_vector

            #print "%10.7f %10.7f %10.7f"%H1.elems

            #print "PRT","%18s"%str((hkl+difference).elems),"idx",omegaidx, "%7.2f"%(
            #  180.*(omegas[omegaidx] + dangles[omegaidx][n]*0.1)/math.pi),

            rot_mat = rotation_vector.axis_and_angle_as_r3_rotation_matrix(
              omegas[omegaidx] + dangles[omegaidx][n]*0.1)
            HP = (rot_mat * matrix.sqr(Amat))*matrix.col(hkl+difference) + beam_vector

            #print "%10.7f %10.7f %10.7f"%HP.elems

            assert approx_equal(H1,HP,eps=1E-5)

            #print


if __name__ == '__main__':

    wavelength = 1.2
    resolution = 3.0
    Amat = (0.0038039968697463817, 0.004498689311366309, 0.0044043429203887785,
      -0.00477859183801569, 0.006594300357213904, -0.002402759536958918,
      -0.01012056453488894, -0.0014226943325514182, 0.002789954423701981)
    orient = crystal_orientation(Amat,True)
    calc = partial_spot_position_partial_H(
      limiting_resolution = resolution,
      orientation = orient.reciprocal_matrix(),
      wavelength = wavelength,
      axial_direction = (0.,1.,0.)
    )
    rotation_scattering(calc,matrix.col((0.,1.,0.)),Amat,wavelength)
    test_finite_differences(calc,matrix.col((0.,1.,0.)),Amat,wavelength)

    print "OK"
