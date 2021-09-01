from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
import math
from libtbx.test_utils import approx_equal
from six.moves import range
from scipy import special
import numpy as np
from simtbx.nanoBragg.tst_gaussian_mosaicity import plotter, plotter2, check_distributions

"""Nick Sauter 6/21/2020
The purpose of function run_uniform() is to provide UMAT and d_UMAT_d_eta with the following
1. Large (scalable) ensemble of UMATs for use by nanoBragg.
2. UMATs are parameterized by eta, the half-angle of mosaic rotation (degrees) assuming Gaussian distribution
3. eta is the standard deviation of the Gaussian
4. The UMATs and d_UMAT_d_eta are completely determined by eta.
5. The UMATs are implemented using the axis & angle approach.
6. Unit axes are drawn randomly from the unit sphere by scitbx.
7. The angle theta is applied, as well as -theta; therefore the number of UMATs is always even.
8. Theta is uniformly chosen over the CDF of the Normal distribution.
9. Therefore the UMAT derivative with respect to eta requires the inverse erf function, supplied by SciPy.
"""

MOSAIC_SPREAD = 2.0 # top hat half width rotation in degrees

def run_uniform(eta_angle, sample_size=20000, verbose=True):

  UMAT = flex.mat3_double()
  d_UMAT_d_eta = flex.mat3_double()
  # the axis is sampled randomly on a sphere.
  # Alternately it could have been calculated on a regular hemispheric grid
  # but regular grid is not needed; the only goal is to have the ability to compute gradient
  mersenne_twister = flex.mersenne_twister(seed=0) # set seed, get reproducible results

  # for each axis, sample both the + and - angle
  assert sample_size%2==0
  # the angle is sampled uniformly from its distribution

  mosaic_rotation0 = np.array(range(sample_size//2))
  mosaic_rotation1 = special.erfinv(mosaic_rotation0/(sample_size//2))
  d_theta_d_eta = math.sqrt(2.0)*flex.double(mosaic_rotation1)
  mosaic_rotation = (math.pi/180.)*eta_angle*d_theta_d_eta

  for m, d in zip(mosaic_rotation, d_theta_d_eta):
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )
    UMAT.append( site.axis_and_angle_as_r3_rotation_matrix(-m,deg=False) )
    d_umat = site.axis_and_angle_as_r3_derivative_wrt_angle(m, deg=False)
    d_UMAT_d_eta.append( (math.pi/180.) * d * d_umat )
    d_umat = site.axis_and_angle_as_r3_derivative_wrt_angle(-m, deg=False)
    d_UMAT_d_eta.append( (math.pi/180.) * -d * d_umat )

  #sanity check on the gaussian distribution
  if verbose:
    nm_angles = check_distributions.get_angular_rotation(UMAT)
    nm_rms_angle = math.sqrt(flex.mean(nm_angles*nm_angles))
    print("Normal rms angle is ",nm_rms_angle)

  return UMAT, d_UMAT_d_eta

def check_finite(mat1, mat2, dmat1, eps):
  # confirms d_Umat_d_eta by finite difference
  for im in range(len(mat1)):
    m1 = sqr(mat1[im])
    m2 = sqr(mat2[im])
    f_diff = (m2 - m1)/eps
    deriv = sqr(dmat1[im])
    assert approx_equal(f_diff.elems, deriv.elems)

def tst_all(make_plots):
  eps =0.00001 # degrees
  UM , d_UM = run_uniform(eta_angle=MOSAIC_SPREAD)
  UMp , d_UMp = run_uniform(eta_angle=MOSAIC_SPREAD+eps)
  UMm , d_UMm = run_uniform(eta_angle=MOSAIC_SPREAD-eps)
  check_finite (UM, UMp, d_UM, eps)
  check_finite (UM, UMm, d_UM, -eps)

  UM_2 , d_UM2 = run_uniform(eta_angle=MOSAIC_SPREAD/2.)
  if make_plots:
    P = plotter(UM,UM_2)
  Q = plotter2(UM,UM_2,make_plots) # suggested by Holton:
  """apply all the UMATs to a particular starting unit vector and take the rms of the resulting end points.
    You should not see a difference between starting with a vector with x,y,z = 0,0,1 vs
    x,y,z = 0.57735,0.57735,0.57735.  But if the implementation is wrong and the UMATs are being made by
    generating Gaussian-random rotations about the three principle axes, then the variance of
    0,0,1 will be significantly smaller than that of 0.57735,0.57735,0.57735."""

if __name__=="__main__":
  import sys
  make_plots = "--plot" in sys.argv
  tst_all(make_plots)
  print("OK")
