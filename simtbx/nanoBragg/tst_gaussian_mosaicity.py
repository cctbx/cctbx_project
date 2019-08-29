from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
import math
from libtbx.test_utils import approx_equal
from six.moves import range

class plotter:
  def __init__(self, tophat, normal):
    # take U-mats from two different distributions, apply them to unit vectors, and plot
    from matplotlib import pyplot as plt
    fig, axes = plt.subplots(2, 3,figsize=(12,7))

    # columns plot the transformation of x, y, and z unit vectors
    rows = [tophat,normal]
    for irow,dist in enumerate(rows):
      iaxes = axes[irow];
      for icol, permutation in enumerate([(0,1,2), (2,0,1), (1,2,0)]):
        axis = iaxes[icol]
        unit = self.permute_vector(vector=(1,0,0), perm = permutation)
        perm_y = self.permute_vector(vector=(0,1,0), perm = permutation).index(1)
        perm_z = self.permute_vector(vector=(0,0,1), perm = permutation).index(1)
        a2 = flex.double(); a3 = flex.double()
        for u in dist:
          U = sqr(u)
          newvec = U * unit
          #print newvec.elems
          a2.append(newvec[perm_y]); a3.append(newvec[perm_z])
        axis.plot (a2,a3,'r,')
        axis.set_aspect("equal")
        axis.set_title("Transformation of unit vector %s"%(str(unit)))
        axis.set_xlim(-0.05,0.05)
        axis.set_ylim(-0.05,0.05)

    plt.show()

  def permute_vector(self, vector, perm):
    return (vector[perm[0]], vector[perm[1]], vector[perm[2]],)

class plotter2:  # compare the transformation of 001 with that of .57735,.57735,.57735
  def __init__(self, tophat, normal, plot):
    # take U-mats from two different distributions, apply them to unit vectors, and plot
    if plot:
      from matplotlib import pyplot as plt
      fig, axes = plt.subplots(2, 2,figsize=(8,7))
    else:
      axes = ((1,2),(3,4)) #dummy

    # columns plot the transformation of x, y, and z unit vectors
    rows = [tophat,normal]
    differences = []
    for irow,dist in enumerate(rows):
      iaxes = axes[irow];
      cube_diag = math.sqrt(1./3) # 0.57735
      for icol, RLP in enumerate([(0,0,1), (cube_diag, cube_diag, cube_diag)]):
        RLP = col(RLP)
        print("(%7.5f %7.5f %7.5f)"%(RLP.elems),"Vector length:%8.6f"%(RLP.length()), end=' ')
        axis = iaxes[icol]
        unit = RLP.normalize()
        seed = col((1,0,0))
        perm2 = unit.cross(seed)
        perm3 = unit.cross(perm2)
        a2 = flex.double(); a3 = flex.double()
        difference_vectors = flex.vec3_double()
        for u in dist:
          U = sqr(u)
          newvec = U * RLP
          difference_vectors.append( newvec-RLP )
          a2.append(newvec.dot(perm2)); a3.append(newvec.dot(perm3))
        rms = math.sqrt( flex.mean ( difference_vectors.dot(difference_vectors) ) )
        print("The rms difference is", rms)
        differences.append(rms)
        if plot:
          axis.plot (a2,a3,'r,')
          axis.set_aspect("equal")
          axis.set_title("Transformation of vector %s"%(str(RLP.elems)))
          axis.set_xlim(-0.05,0.05)
          axis.set_ylim(-0.05,0.05)

      assert approx_equal(differences[0],differences[1],eps=1e-04), \
      "RMS mosaic distribution for axis vector and diagonal vector should be similar, as proposed by J Holton"

    if plot: plt.show()

class check_distributions:
  @staticmethod
  def get_angular_rotation(dist):
    angle_deg = flex.double()
    for umat in dist:
      angle,axis = sqr(
        umat).r3_rotation_matrix_as_unit_quaternion().unit_quaternion_as_axis_and_angle(deg=True)
      angle_deg.append(angle)
    return angle_deg

MOSAIC_SPREAD = 2.0 # top hat half width rotation in degrees
SAMPLE_SIZE = 90000

def run_sim2smv(fileout):
  SIM = nanoBragg(detpixels_slowfast=(1000,1000),pixel_size_mm=0.1,Ncells_abc=(5,5,5),verbose=0)
  SIM.mosaic_spread_deg = MOSAIC_SPREAD # apparently this is half width
  SIM.mosaic_domains = SAMPLE_SIZE
  SIM.distance_mm=100 # this triggers the generation of mosaic distribution
  UMAT_th = SIM.get_mosaic_blocks() # extract top-hat distributed U-mats

  #sanity checks on the top hat distribution
  th_angles = check_distributions.get_angular_rotation(UMAT_th)
  max_angle = flex.max(th_angles)
  assert max_angle <= MOSAIC_SPREAD + 0.0000001 # need to allow a small epsilon
  assert max_angle > 0.99 * MOSAIC_SPREAD # insist that max angle is near the limit we gave it
  rms_angle = math.sqrt(flex.mean(th_angles*th_angles))
  print(rms_angle)
  assert rms_angle < MOSAIC_SPREAD

  import scitbx # compute an array of normally-distributed U-mats into the simulator
  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0) # set seed, get reproducible results
  scitbx.random.set_random_seed(4321) # set seed, get reproducibe results
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=rms_angle*math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(SAMPLE_SIZE)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )

  # set the normally-distributed U-mats into the simulator:
  SIM.set_mosaic_blocks(UMAT_nm)

  # get them back
  new_UMAT_nm = SIM.get_mosaic_blocks()

  #double check that simulator does not alter the UMATs
  for iumat in range(0,SAMPLE_SIZE,100):




    assert UMAT_nm[iumat]==new_UMAT_nm[iumat]

  #sanity check on the gaussian distribution
  nm_angles = check_distributions.get_angular_rotation(new_UMAT_nm)
  nm_rms_angle = math.sqrt(flex.mean(nm_angles*nm_angles))
  print(nm_rms_angle)
  assert approx_equal(rms_angle,nm_rms_angle,eps=1e-03), \
  "The top hat and gaussian models should have similar standard deviations"

  return UMAT_th, new_UMAT_nm

def tst_all(make_plots):
  fileout = "dummy_file.cbf"
  UM_th, UM_nm = run_sim2smv(fileout)
  if make_plots:
    P = plotter(UM_th,UM_nm)
  Q = plotter2(UM_th,UM_nm,make_plots) # suggested by Holton:
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
