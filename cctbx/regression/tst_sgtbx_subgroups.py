from __future__ import absolute_import, division, print_function
from cctbx.sgtbx import subgroups
from cctbx import sgtbx
from six.moves import range

class subgroup_stats(object):

  def __init__(self, parent_group_info):
    self.subgroups = subgroups.subgroups(
      parent_group_info).groups_parent_setting()
    self.n_non_centric = 0
    self.n_chiral = 0
    for subgroup in self.subgroups:
      if (not subgroup.is_centric()): self.n_non_centric += 1
      if (subgroup.is_chiral()): self.n_chiral += 1

def run():
  for space_group_symbol in ("P-1",
                             "P2/m",
                             "C2/m",
                             "Pmmm",
                             "Cmmm",
                             "Fmmm",
                             "Immm",
                             "P4/mmm",
                             "I4/mmm",
                             "R-3m",
                             "P6/mmm",
                             "Pm-3m",
                             "Im-3m",
                             "Fm-3m"):
    centric_info = sgtbx.space_group_info(space_group_symbol)
    non_centric = sgtbx.space_group()
    for i_ltr in range(centric_info.group().n_ltr()):
      for i_smx in range(centric_info.group().n_smx()):
        s = centric_info.group()(i_ltr,0,i_smx)
        non_centric.expand_smx(s)
    assert non_centric.f_inv() == 1
    assert non_centric.order_z() * 2 == centric_info.group().order_z()
    non_centric_info = sgtbx.space_group_info(group=non_centric)
    centric_stats = subgroup_stats(centric_info)
    non_centric_stats = subgroup_stats(non_centric_info)
    assert len(centric_stats.subgroups) >= 2*len(non_centric_stats.subgroups)
    assert centric_stats.n_non_centric >= non_centric_stats.n_non_centric
    assert centric_stats.n_chiral == non_centric_stats.n_chiral
    assert non_centric_stats.n_non_centric == len(non_centric_stats.subgroups)
    assert non_centric_stats.n_non_centric == non_centric_stats.n_chiral
  print("OK")

if (__name__ == "__main__"):
  run()
