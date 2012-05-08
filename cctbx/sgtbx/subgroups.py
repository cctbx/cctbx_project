from cctbx import sgtbx

def anomalous_reflection_intensity_primitive_cell(space_group):
  assert space_group.n_ltr() == 1
  assert not space_group.is_centric()
  subgroup = sgtbx.space_group(space_group)
  subgroup.make_tidy()
  result = [subgroup]
  for s1 in space_group:
    assert s1.t().num() == (0,0,0)
    for s2 in space_group:
      subgroup = sgtbx.space_group()
      subgroup.expand_smx(s1)
      subgroup.expand_smx(s2)
      subgroup.make_tidy()
      if (not subgroup in result):
        result.append(subgroup)
  return result

class subgroups(object):

  def __init__(self, parent_group_info):
    self._p_groups = []
    self.z2p_op = parent_group_info.group().z2p_op()
    p_parent_group_info = parent_group_info.change_basis(self.z2p_op)
    p_parent_group = p_parent_group_info.group()
    assert p_parent_group.order_p() == p_parent_group.order_z()
    p_parent_group.make_tidy()
    for i_smx in xrange(p_parent_group.order_p()):
      group_i = sgtbx.space_group()
      group_i.expand_smx(p_parent_group(i_smx))
      for j_smx in xrange(i_smx,p_parent_group.order_p()):
        subgroup = sgtbx.space_group(group_i)
        subgroup.expand_smx(p_parent_group(j_smx))
        subgroup.make_tidy()
        self._add(subgroup)

  def _add(self, group):
    for g in self._p_groups:
      if (g == group): return 0
    self._p_groups.append(group)
    return 1

  def groups_primitive_setting(self):
    return self._p_groups

  def groups_parent_setting(self):
    result = []
    p2z_op = self.z2p_op.inverse()
    for g in self._p_groups:
      result.append(g.change_basis(p2z_op))
    return result

def show(parent_group_info):
  parent_group_info.show_summary()
  subgrs = subgroups(parent_group_info).groups_parent_setting()
  print "number of subgroups:", len(subgrs)
  for subgroup in subgrs:
    subgroup_info = sgtbx.space_group_info(group=subgroup)
    subgroup_info.show_summary()
  print

if (__name__ == "__main__"):
  raise RuntimeError("Please use the cctbx.subgroups command.")
