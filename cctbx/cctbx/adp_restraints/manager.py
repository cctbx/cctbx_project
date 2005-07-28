from cctbx import crystal
from cctbx.array_family import flex
import math
import os

class iso:
  def __init__(self, xray_structure,
                     geometry_restraints_manager,
                     sphere_radius,
                     distance_power,
                     wilson_b,
                     mean_power,
                     normalize,
                     collect=False):
     assert normalize in (True, False)
     grm = geometry_restraints_manager
     assert grm.plain_pair_sym_table is not None
     assert grm.plain_pairs_radius is not None
     assert sphere_radius <= grm.plain_pairs_radius
     uc = xray_structure.unit_cell()
     sites = xray_structure.sites_frac()
     u_isos = xray_structure.extract_u_iso_or_u_equiv()
     self._target = 0.0
     self._gradients = flex.double(sites.size(),0.0)
     self.counter = 0
     self.obj = crystal.adp_iso_restraint_helper(
                      pair_sym_table           = grm.plain_pair_sym_table,
                      orthogonalization_matrix = uc.orthogonalization_matrix(),
                      sites_frac               = sites,
                      u_isos                   = u_isos,
                      sphere_radius            = sphere_radius,
                      distance_power           = distance_power,
                      mean_power               = mean_power,
                      normalize                = normalize,
                      collect                  = collect)
     self._target = self.obj.target()
     self._gradients = self.obj.derivatives()
     self.counter = self.obj.number_of_members()
     if(wilson_b is not None):
        u_mean = flex.mean(u_isos)
        u_wilson = wilson_b/(math.pi**2*8)
        term = u_mean - u_wilson
        gw = 2.*term / u_isos.size()
        gwa = flex.double(u_isos.size(), gw)
        w = (math.sqrt(flex.sum(flex.pow2(self._gradients))) / \
             math.sqrt(flex.sum(flex.pow2(gwa))))
        self._gradients = self._gradients * w + gwa
        self._target = self._target * w + term**2

  def target(self):
    return self._target

  def gradients(self):
    return self._gradients
