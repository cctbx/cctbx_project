from cctbx import adptbx
from scitbx.python_utils.misc import adopt_init_args

def quality_factor_from_any(d_min=None,
                            grid_resolution_factor=None,
                            quality_factor=None,
                            u_extra=None,
                            b_extra=None):
  assert [quality_factor, u_extra, b_extra].count(None) >= 2
  if (u_extra is not None):
    b_extra = adptbx.u_as_b(u_extra)
  if (b_extra is not None):
    assert [d_min, grid_resolution_factor].count(None) == 0
    assert d_min > 0
    sigma = 1 / (2. * grid_resolution_factor)
    log_quality_factor = b_extra * sigma * (sigma - 1) / (d_min * d_min)
    quality_factor = 10**log_quality_factor
  elif (quality_factor is None):
    quality_factor = 100
  return quality_factor

class from_scatterers_common:

  def __init__(self, manager, xray_structure, miller_set):
    adopt_init_args(self, locals(), hide=0001)
    assert xray_structure is not None and miller_set is not None
    assert xray_structure.unit_cell().is_similar_to(miller_set.unit_cell())
    assert xray_structure.space_group() == miller_set.space_group()
    if (manager is not None):
      assert xray_structure.unit_cell().is_similar_to(manager.unit_cell())
      assert xray_structure.space_group() == manager.space_group()

  def manager(self):
    return self._manager

  def xray_structure(self):
    return self._xray_structure

  def miller_set(self):
    return self._miller_set
