from __future__ import absolute_import, division, print_function
from cctbx.development import random_structure
from cctbx.sgtbx import space_group_info
from cctbx_asymmetric_map_ext import *
from cctbx.array_family import flex
from cctbx import maptbx
from six.moves import range

result = flex.double()

def run_group(symbol):
  group = space_group_info(symbol);
  elements = ('C',)*10
  xrs = random_structure.xray_structure(
    space_group_info = group,
    volume_per_atom = 50.,
    general_positions_only = True,
    elements = elements,
    min_distance = 2.0,
    u_iso=0.1)
  d_min = 2.
  fc = xrs.structure_factors(d_min=d_min).f_calc()
  symmetry_flags = None#maptbx.use_space_group_symmetry
  fft_map = fc.fft_map(symmetry_flags = symmetry_flags, resolution_factor=1./5)
  fft_map.apply_volume_scaling()
  map_data = fft_map.real_map_unpadded()
  #
  r = maptbx.peak_volume_estimate(
    map_data         = map_data,
    sites_cart       = xrs.sites_cart(),
    crystal_symmetry = xrs.crystal_symmetry(),
    cutoff           = 0.3,
    atom_radius      = 1.5)
  result.append(r)

def run():
  for i in range(1,231):
    run_group(i);

if (__name__ == "__main__"):
  run()
  print(result.min_max_mean().as_tuple()) # XXX why they are so different ?
