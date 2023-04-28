from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex

import boost_adaptbx.boost.python as bp
from six.moves import zip
ext = bp.import_ext("cctbx_other_restraints_ext")
from cctbx_other_restraints_ext import *

from cctbx.geometry_restraints import weight_as_sigma
import sys

@bp.inject_into(sump_proxy)
class _():

  def show(self, site_occupancies, site_labels, f, prefix):
    r = sump(site_occupancies, self)
    label = []
    for i, i_seq in enumerate(self.i_seqs):
      if self.labels:
        l = self.labels[i]
      else:
        l = "%s.occu" %(site_labels[i_seq] if site_labels else str(i_seq))
      k = self.coefficients[i]
      if abs(k-1) < 1e-2:
        p = l
      elif abs(k-round(k)) < 1e-2:
        p = "%s*%s" %(round(k), l)
      else:
        p = "%.2f*%s" %(self.coefficients[i], l)
      if label and p[0] != '-':
        p = '+ ' + p
      label.append(p)
    label = ' '.join(label)
    label += " = %.2f" %r.target
    print("%s          %s delta=%.2e    sigma=%.2e   weight=%.2e"\
      %(prefix, label, r.delta, weight_as_sigma(r.weight), r.weight), end=' ', file=f)

@bp.inject_into(shared_sump_proxy)
class _():

  def show(self,
        site_occupancies,
        site_labels=None,
        f=None,
        prefix=""):
    for proxy in self:
      proxy.show(site_occupancies, site_labels, f, prefix)

  def linearise(self,
        site_occupancies,
        param_map,
        linearised_eqns):
    for proxy in self:
      r = sump(site_occupancies, proxy)
      r.linearise(linearised_eqns, param_map, proxy.i_seqs)


