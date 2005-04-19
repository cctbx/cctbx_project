from cctbx.array_family import flex
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import maptbx
from cctbx.maptbx.real_space_refinement import lbfgs
from cctbx import uctbx
from scitbx.python_utils.math_utils import ifloor, iceil
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
from cStringIO import StringIO
import math
import sys

def exercise_real_space_refinement(verbose):
  if (verbose):
    out = sys.stdout
  else:
    out = StringIO()
  crystal_symmetry = crystal.symmetry(
    unit_cell=(10,10,10,90,90,90),
    space_group_symbol="P 1")
  xray_structure = xray.structure(
    crystal_symmetry=crystal_symmetry,
    scatterers=flex.xray_scatterer([
      xray.scatterer(label="C", site=(0,0,0))]))
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry,
    anomalous_flag=False,
    d_min=1)
  f_calc = miller_set.structure_factors_from_scatterers(
    xray_structure=xray_structure).f_calc()
  fft_map = f_calc.fft_map()
  fft_map.apply_sigma_scaling()
  real_map = fft_map.real_map_unpadded()
  #
  minfrac = crystal_symmetry.unit_cell().fractionalize((-5,-5,-5))
  maxfrac = crystal_symmetry.unit_cell().fractionalize((5,5,5))
  gridding_first = [ifloor(n*b) for n,b in zip(fft_map.n_real(), minfrac)]
  gridding_last = [iceil(n*b) for n,b in zip(fft_map.n_real(), maxfrac)]
  data=maptbx.copy(real_map, gridding_first, gridding_last)
  #
  grid_cell=uctbx.unit_cell((10.0/30,10.0/30,10.0/30,90,90,90))
  grid_mat = grid_cell.fractionalization_matrix()
  delta_h = .005
  for i_trial in xrange(100):
    sites_cart = flex.vec3_double((flex.random_double(size=3)-0.5)*1)
    tmp_sites_cart = sites_cart.deep_copy()
    for i in xrange(3):
      ref = lbfgs(
        data=data,
        grid_mat=grid_mat,
        sites_cart=tmp_sites_cart,
        delta_h=delta_h)
      temp = flex.double(ref.sites_cart[0])-flex.double((0,0,0))
      temp = math.sqrt(temp.dot(temp))
      if temp <= 2*delta_h:
        break
      print >> out, "recycling:", ref.sites_cart[0]
      tmp_sites_cart = ref.sites_cart
    for site,sitec in zip(ref.sites_cart,xray_structure.sites_cart()):
      print >> out, i_trial
      print >> out, sitec
      print >> out, sites_cart[0]
      print >> out, site
      temp = flex.double(site)-flex.double(sitec)
      temp = math.sqrt(temp.dot(temp))
      print >> out, temp, delta_h
      assert temp <= delta_h*2
      print >> out

if (__name__ == "__main__"):
  exercise_real_space_refinement("--verbose" in sys.argv[1:])
  print format_cpu_times()
