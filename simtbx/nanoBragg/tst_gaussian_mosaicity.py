from __future__ import division
from scitbx.array_family import flex
from scitbx.matrix import sqr,col
from simtbx.nanoBragg import nanoBragg
import libtbx.load_env # possibly implicit
import math

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

def run_sim2smv(fileout):
  SIM = nanoBragg(detpixels_slowfast=(1000,1000),pixel_size_mm=0.1,Ncells_abc=(5,5,5),verbose=0)
  SIM.mosaic_domains = 10000
  SIM.mosaic_spread_deg = MOSAIC_SPREAD
  SIM.distance_mm=100 # this triggers the generation of mosaic distribution
  UMAT_th = SIM.get_mosaic_blocks() # extract top-hat distributed U-mats

  #sanity checks on the top hat distribution
  th_angles = check_distributions.get_angular_rotation(UMAT_th)
  max_angle = flex.max(th_angles)
  assert max_angle <= MOSAIC_SPREAD + 0.0000001 # need to allow a small epsilon
  assert max_angle > 0.99 * MOSAIC_SPREAD # insist that max angle is near the limit we gave it
  rms_angle = math.sqrt(flex.mean(th_angles*th_angles))
  print rms_angle

  import scitbx # computer an array of normally-distributed U-mats into the simulator
  UMAT_nm = flex.mat3_double()
  mersenne_twister = flex.mersenne_twister(seed=0)
  rand_norm = scitbx.random.normal_distribution(mean=0, sigma=rms_angle*math.pi/180.)
  g = scitbx.random.variate(rand_norm)
  mosaic_rotation = g(10000)
  for m in mosaic_rotation:
    site = col(mersenne_twister.random_double_point_on_sphere())
    UMAT_nm.append( site.axis_and_angle_as_r3_rotation_matrix(m,deg=False) )

  # set the normally-distributed U-mats into the simulator:
  SIM.set_mosaic_blocks(UMAT_nm)

  # get them back
  new_UMAT_nm = SIM.get_mosaic_blocks()

  #sanity check on the gaussian distribution
  nm_angles = check_distributions.get_angular_rotation(new_UMAT_nm)
  nm_rms_angle = math.sqrt(flex.mean(nm_angles*nm_angles))
  print nm_rms_angle

  return UMAT_th, new_UMAT_nm

def tst_all(make_plots):
  fileout = "dummy_file.cbf"
  UM_th, UM_nm = run_sim2smv(fileout)
  if make_plots:
    P = plotter(UM_th,UM_nm)

if __name__=="__main__":
  import sys
  make_plots = "--plot" in sys.argv
  tst_all(make_plots)
  print "OK"
