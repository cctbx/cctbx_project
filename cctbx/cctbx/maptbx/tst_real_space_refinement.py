from cctbx.array_family import flex
from cctbx import xray
from cctbx import miller
from cctbx import crystal
from cctbx import maptbx
from cctbx.maptbx.real_space_refinement import lbfgs
from cctbx import uctbx
from scitbx.python_utils.math_utils import ifloor, iceil
from libtbx.test_utils import approx_equal

def exercise_real_space_refinement():
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
  for i_trial in xrange(100):
    sites_cart = flex.vec3_double((flex.random_double(size=3)-0.5)*1)
    ref = lbfgs(
      data=data,
      grid_mat=grid_mat,
      sites_cart=sites_cart.deep_copy())
    for site,sitec in zip(ref.sites_cart,xray_structure.sites_cart()):
      print sitec
      print sites_cart[0]
      print site
      assert approx_equal(site,sitec)
      print

if (__name__ == "__main__"):
  exercise_real_space_refinement()
  print "OK"
