from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
import os

from iotbx_boost import scalepack_ext

class no_merge_original_index_file:

  # scalepack manual, edition 5, page 132
  # no merge
  # original index
  # the output will also contain the original (not unique) hkl for each
  # reflection. This is designed for MAD/local scaling work. The
  # original index modifier only works with the default INCLUDE
  # NO PARTIALS. The output will consist of the original hkl, unique
  # hkl, batch number, a flag (0 = centric, 1 = I+, 2 = I-), another flag
  # (0 = hkl reflecting above the spindle, 1 = hkl reflecting below the
  # spindle), the asymmetric unit of the reflection, I (scaled, Lorentz
  # and Polarization corrected), and the s of I. The format is
  # (6i4, i6, 2i2, i3, 2f8.1).
  # 
  #    orig. hkl   uniq. hkl    b# c s  a       I    sigI
  #    0   0   3   0   0   3    14 0 0  1    -1.8     1.3

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
    assert self.space_group.order_p() == n_sym_ops_from_file
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

  def merge_equivalents(self, crystal_symmetry=None):
    if (crystal_symmetry is None):
      crystal_symmetry = crystal.symmetry(space_group=self.space_group)
    return miller.array(
      miller_set=miller.set(
        crystal_symmetry=crystal_symmetry,
        indices=self.original_indices,
        anomalous_flag=0001),
      data=self.i_obs,
      sigmas=self.sigmas).merge_equivalents()

def quick_test(file_name):
  from scitbx.python_utils.misc import user_plus_sys_time
  t = user_plus_sys_time()
  s = no_merge_original_index_file(file_name)
  print "Time read:", t.delta()
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
  m = s.merge_equivalents()
  m.array().show_summary()
  print "min redundancies:", flex.min(m.redundancies())
  print "max redundancies:", flex.max(m.redundancies())
  print "mean redundancies:", flex.mean(m.redundancies().as_double())
  print

def run():
  file_names = """
/net/redbelly/lbnl1/share/structure-lib/p9/data/infl.sca
/net/redbelly/lbnl1/share/structure-lib/p9/data/peak.sca
/net/redbelly/lbnl1/share/structure-lib/p9/data/high.sca
/net/redbelly/lbnl1/share/structure-lib/rh-dehalogenase/data/auki_rd_1.sca
/net/redbelly/lbnl1/share/structure-lib/rh-dehalogenase/data/hgi2_rd_1.sca
/net/redbelly/lbnl1/share/structure-lib/rh-dehalogenase/data/hgki_rd_1.sca
/net/redbelly/lbnl1/share/structure-lib/rh-dehalogenase/data/ndac_rd_1.sca
/net/redbelly/lbnl1/share/structure-lib/rh-dehalogenase/data/rt_rd_1.sca
/net/redbelly/lbnl1/share/structure-lib/rh-dehalogenase/data/smac_1.sca
/net/redbelly/lbnl1/share/structure-lib/vmp/data/infl.sca
/net/redbelly/lbnl1/share/structure-lib/vmp/data/peak.sca
/net/redbelly/lbnl1/share/structure-lib/vmp/data/high.sca
/net/redbelly/scratch1/rwgk/bnl_2003/karen/shelxd/p123-unmerged.sca
/net/redbelly/scratch1/rwgk/bnl_2003/karen/shelxd/pk1-unmerged.sca
/net/redbelly/scratch1/rwgk/bnl_2003/karen/shelxd/pk12-unmerged.sca
/net/redbelly/scratch1/rwgk/bnl_2003/karen/shelxd/pk1234-unmerged.sca""".split()
  for file_name in file_names:
    print "File name:", file_name
    quick_test(file_name)

if (__name__ == "__main__"):
  run()
