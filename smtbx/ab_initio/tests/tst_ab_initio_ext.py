from __future__ import absolute_import, division, print_function
from cctbx import maptbx
from cctbx.development import random_structure
from cctbx import miller
from cctbx import sgtbx

from scitbx.array_family import flex

from  libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times

from smtbx import ab_initio

def run():
  structure = random_structure.xray_structure(
    sgtbx.space_group_info("P21/c"),
    elements=["Si"]*10,
    volume_per_atom=18.6,
    min_distance=1.2,
    general_positions_only=False)
  miller_set_f_obs = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=True,
    d_min=0.8)
  f_obs = miller_set_f_obs.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc()
  fft_map = f_obs.fft_map(symmetry_flags=maptbx.use_space_group_symmetry)

  padded = fft_map.real_map()
  unpadded = fft_map.real_map_unpadded() # copy
  unpadded_1d = unpadded.as_1d() # 1D view => in-place
  mmm = flex.min_max_mean_double(unpadded_1d)

  for delta in ((mmm.min + mmm.mean)/2, mmm.mean, (mmm.mean + mmm.max)/2):
    # in-place charge flipping
    ab_initio.ext.flip_charges_in_place(padded, delta)

    # same but on an unpadded copy using the flex tools
    flipped_selection = unpadded_1d < delta
    flipped = unpadded_1d.select(flipped_selection)
    flipped *= -1
    unpadded_1d.set_selected(flipped_selection, flipped)

    assert approx_equal(padded, unpadded, 1e-15)

  print(format_cpu_times())


if __name__ == '__main__':
  run()
