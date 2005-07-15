from mmtbx import monomer_library
import mmtbx.monomer_library.server
import mmtbx.monomer_library.pdb_interpretation
import sys, math
from libtbx.test_utils import approx_equal
from iotbx import pdb
from cctbx.array_family import flex
from mmtbx.monomer_library import pdb_interpretation
from cctbx.geometry_restraints.lbfgs import lbfgs as geometry_restraints_lbfgs
import scitbx.lbfgs
import libtbx.load_env
import os

class iso:
  def __init__(self, xray_structure,
                     geometry_restraints_manager,
                     sphere_radius,
                     power_factor,
                     wilson_b=None,
                     p=0,
                     normalize=True):
     assert p == 1 or p == 0
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

     for i_seq,pair_sym_dict in enumerate(grm.plain_pair_sym_table):
         for j_seq,sym_ops in pair_sym_dict.items():
             for sym_op in sym_ops:
                 d = uc.distance(sites[i_seq], sym_op*sites[j_seq])
                 if(d <= sphere_radius and d > 0.0):
                    weight = 1/d**power_factor
                    u1 = u_isos[i_seq]
                    u2 = u_isos[j_seq]
                    term = u1 - u2
                    if(p == 0):
                       self._gradients[i_seq] += weight * term * 2
                       self._gradients[j_seq] -= weight * term * 2
                       self._target += weight * term**2
                    else:
                       sum = u1 + u2
                       g1 = (u1**2-3*u2**2+2*u1*u2) / sum**2
                       g2 = (u2**2-3*u1**2+2*u1*u2) / sum**2
                       self._gradients[i_seq] += weight * g1
                       self._gradients[j_seq] += weight * g2
                       self._target += weight * term**2 / sum
                    self.counter += 1
     if(normalize):
        self._target /= self.counter
        self._gradients = self._gradients / self.counter
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
