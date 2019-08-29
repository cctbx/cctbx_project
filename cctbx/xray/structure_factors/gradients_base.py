from __future__ import absolute_import, division, print_function
from cctbx.xray.structure_factors.manager import managed_calculation_base

class gradients_base(managed_calculation_base):

  def d_target_d_f_calc(self):
    return self._d_target_d_f_calc

  def packed(self):
    return self._results.packed()

  def check_size(self, array):
    assert self.packed().size() == 0
    assert array.size() == self.xray_structure().scatterers().size()
    return array

  def d_target_d_u_iso(self):
    return self.check_size(self._results.d_target_d_u_iso())

  def d_target_d_occupancy(self):
    return self.check_size(self._results.d_target_d_occupancy())

  def d_target_d_fp(self):
    return self.check_size(self._results.d_target_d_fp())

  def d_target_d_fdp(self):
    return self.check_size(self._results.d_target_d_fdp())
