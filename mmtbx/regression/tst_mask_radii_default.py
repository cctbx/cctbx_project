from __future__ import absolute_import, division, print_function
"""
Default bulk-solvent mask radii. The values were derived (2026-11) on 55
low-resolution and 10 high-resolution refinements with and without riding H:
0.95/0.70 is non-inferior to the historical 1.1/0.9 without hydrogens and
better with them (see the comment at mask_master_params in mmtbx.masks).
"""
import iotbx.phil
import mmtbx.masks
import mmtbx.f_model
from cctbx import sgtbx
from cctbx.development import random_structure
from libtbx.test_utils import approx_equal

def test_phil_defaults():
  p = mmtbx.masks.mask_master_params.extract()
  assert approx_equal(p.solvent_radius, 0.95)
  assert approx_equal(p.shrink_truncation_radius, 0.70)
  # Plain numbers: the defaults do not depend on the model.
  assert isinstance(p.solvent_radius, float)
  assert isinstance(p.shrink_truncation_radius, float)
  # Round trip through the PHIL text.
  s = mmtbx.masks.mask_master_params.format(python_object=p).as_str()
  assert "solvent_radius = 0.95" in s, s
  assert "shrink_truncation_radius = 0.7" in s, s
  # Explicit values are kept.
  user = iotbx.phil.parse("solvent_radius=1.1\nshrink_truncation_radius=0.9")
  q = mmtbx.masks.mask_master_params.fetch(source=user).extract()
  assert approx_equal((q.solvent_radius, q.shrink_truncation_radius),
    (1.1, 0.9))

def test_fmodel_uses_defaults():
  xrs = random_structure.xray_structure(
    space_group_info       = sgtbx.space_group_info("P 21 21 21"),
    elements               = ("C","N","O")*8 + ("H",)*12,
    volume_per_atom        = 50,
    min_distance           = 1.0,
    general_positions_only = True,
    random_u_iso           = True)
  f_obs = abs(xrs.structure_factors(d_min=2.5).f_calc())
  flags = f_obs.generate_r_free_flags(fraction=0.1)
  for data_type in ["xray", "neutron"]:
    fm = mmtbx.f_model.manager(f_obs=f_obs, r_free_flags=flags,
      xray_structure=xrs, data_type=data_type)
    assert approx_equal(fm.mask_params.solvent_radius, 0.95)
    assert approx_equal(fm.mask_params.shrink_truncation_radius, 0.70)
    assert approx_equal(fm.mask_manager.mask_params.solvent_radius, 0.95)
    fm.update_all_scales(remove_outliers=False, refine_hd_scattering=False)
    info = fm.info()
    assert approx_equal(info.mask_solvent_radius, 0.95)
    assert approx_equal(info.mask_shrink_radius, 0.70)

if (__name__ == "__main__"):
  test_phil_defaults()
  test_fmodel_uses_defaults()
  print("OK")
