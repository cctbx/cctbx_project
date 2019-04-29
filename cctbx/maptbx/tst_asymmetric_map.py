from __future__ import division, print_function
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
import boost.python
ext = boost.python.import_ext("cctbx_asymmetric_map_ext")
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import maptbx

def run_group(symbol):
  group = space_group_info(symbol);
  print("\n==")
  elements = ('C', 'N', 'O', 'H')*11
  struc = random_structure.xray_structure(
    space_group_info = group,
    volume_per_atom = 25.,
    general_positions_only = False,
    elements = elements,
    min_distance = 1.0)
  struc.show_summary()
  d_min = 2.
  fc = struc.structure_factors(d_min=d_min).f_calc()
  symmetry_flags = maptbx.use_space_group_symmetry
  fftmap = fc.fft_map(symmetry_flags = symmetry_flags)
  grid_size = fftmap.real_map().accessor().focus()

  # grid_tags = maptbx.grid_tags(grid_size)
  # grid_tags.build(fftmap.space_group_info().type(), fftmap.symmetry_flags())
  # print "grid_tags:\n", dir(grid_tags)
  # grid_tags.verify(real_map)

  print("FFT grid_size = ", grid_size)
  amap = asymmetric_map(struc.space_group().type(), fftmap.real_map())
  afc = amap.structure_factors(fc.indices())
  afftmap = amap.map_for_fft()
  print("whole cell map size: ", afftmap.accessor().focus())
  adata = amap.data()
  acc = adata.accessor()
  print("Asu map size: ", acc.origin(), " ", acc.last(), " ", acc.focus(), \
      " ", acc.all())
  df = flex.abs(afc - fc.data())
  r1 = flex.sum(df) / flex.sum(flex.abs(fc.data()))
  print("R1: ", r1)
  assert r1 < 1.e-5
  # just to prove to myself that I can shift origin to 000 and then reshape back
  adata = adata.shift_origin()
  adata.reshape(acc)
  #
  adata2 = adata.deep_copy()*2.
  amap2 = asymmetric_map(struc.space_group().type(), adata2, grid_size)
  afc2 = amap2.structure_factors(fc.indices())
  df2 = flex.abs(afc2*.5-fc.data())
  r12 = flex.sum(df2) / flex.sum(flex.abs(fc.data()))
  print("R1 again: ", r12)
  assert r12 < 1.e-5

  p1_map = amap.symmetry_expanded_map()
  assert p1_map.accessor().focus() == grid_size

  rel_tol = 1.e-6
  n = 0
  mean_rel_dif = 0.
  for (m1,m2) in zip(fftmap.real_map(),p1_map):
    dif = abs(m1-m2)
    av = 0.5*(abs(m1)+abs(m2))
    assert dif <= rel_tol * av, "%f not <= %f * %f" % (dif, rel_tol, av)
    if av != 0:
      mean_rel_dif = mean_rel_dif + dif/av
      n = n + 1
  mean_rel_dif = mean_rel_dif / n
  print("mean rel err: ", mean_rel_dif)
  assert mean_rel_dif < 1.e-6

def run():
  for i in xrange(1,231):
    run_group(i);

if (__name__ == "__main__"):
  run()
