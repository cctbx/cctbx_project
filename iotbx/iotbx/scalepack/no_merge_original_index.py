from cctbx import sgtbx
from cctbx.array_family import flex
import os

from iotbx_boost import scalepack_ext

class no_merge_original_index_file:

  def __init__(self, file_name):
    self.file_name = os.path.normpath(file_name)
    f = open(file_name)
    line = f.readline()
    n_sym_ops_from_file = int(line[:5].strip())
    assert line[5] == " "
    self.space_group_symbol = line[6:].strip()
    self.space_group = sgtbx.space_group()
    for i in xrange(n_sym_ops_from_file):
      line = f.readline()
      r = sgtbx.rot_mx([int(line[j*3:(j+1)*3]) for j in xrange(9)], 1)
      line = f.readline()
      t = sgtbx.tr_vec([int(line[j*3:(j+1)*3]) for j in xrange(3)], 12)
      self.space_group.expand_smx(sgtbx.rt_mx(r, t))
    assert self.space_group.order_z() == n_sym_ops_from_file
    f.close()
    all_arrays = scalepack_ext.no_merge_original_index_arrays(
      file_name, n_sym_ops_from_file*2+1)
    self.original_indices = all_arrays.original_indices()
    self.unique_indices = all_arrays.unique_indices()
    self.batch_numbers = all_arrays.batch_numbers()
    self.centric_tags = all_arrays.centric_tags()
    self.spindle_flags = all_arrays.spindle_flags()
    self.asymmetric_unit_indices = all_arrays.asymmetric_unit_indices()
    self.i_obs = all_arrays.i_obs()
    self.sigmas = all_arrays.sigmas()

from scitbx.python_utils.misc import user_plus_sys_time
t = user_plus_sys_time()
#s = no_merge_original_index_file("head_tail.sca")
s = no_merge_original_index_file("pk1-unmerged.sca")
#s = no_merge_original_index_file("tt.sca")
print t.delta()
print tuple(s.original_indices[:3])
print tuple(s.unique_indices[:3])
print tuple(s.batch_numbers[:3])
print tuple(s.centric_tags[:3])
print tuple(s.spindle_flags[:3])
print tuple(s.asymmetric_unit_indices[:3])
print tuple(s.i_obs[:3])
print tuple(s.sigmas[:3])
print tuple(s.original_indices[-3:])
print tuple(s.unique_indices[-3:])
print tuple(s.batch_numbers[-3:])
print tuple(s.centric_tags[-3:])
print tuple(s.spindle_flags[-3:])
print tuple(s.asymmetric_unit_indices[-3:])
print tuple(s.i_obs[-3:])
print tuple(s.sigmas[-3:])
