from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.str_utils import show_string
import sys, os

import boost.python
scalepack_ext = boost.python.import_ext("iotbx_scalepack_ext")

class reader(object):

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

  def __init__(self, file_name, header_only=False):
    self.file_name = os.path.normpath(file_name)
    f = open(file_name)
    line = f.readline()
    assert line[5] == " "
    n_sym_ops_from_file = int(line[:5].strip())
    assert n_sym_ops_from_file > 0
    self.space_group_symbol = line[6:].strip()
    self.space_group_from_ops = sgtbx.space_group()
    for i in xrange(n_sym_ops_from_file):
      line = f.readline().rstrip()
      assert len(line) == 27
      r = sgtbx.rot_mx([int(line[j*3:(j+1)*3]) for j in xrange(9)], 1)
      line = f.readline().rstrip()
      assert len(line) == 9
      t = sgtbx.tr_vec([int(line[j*3:(j+1)*3]) for j in xrange(3)], 12)
      self.space_group_from_ops.expand_smx(sgtbx.rt_mx(r, t))
    f.close()
    if (header_only):
      self.original_indices = None
      return
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

  def space_group_info(self):
    return combine_symops_and_symbol(
      space_group_from_ops=self.space_group_from_ops,
      space_group_symbol=self.space_group_symbol)

  def show_summary(self, out=None, prefix=""):
    if (out is None): out = sys.stdout
    print >> out, prefix + "File name:", show_string(self.file_name)
    print >> out, prefix + "Space group symbol:", \
      show_string(self.space_group_symbol)
    try: space_group_info = self.space_group_info()
    except KeyboardInterrupt: raise
    except: pass
    else:
      space_group_info.show_summary(
        f=out, prefix=prefix+"Space group from operations: ")
    if (self.original_indices is not None):
      print >> out, prefix + "Number of original indices:", \
        self.original_indices.size()

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=None,
      space_group_info=self.space_group_info())

  def unmerged_miller_set(self, crystal_symmetry=None, force_symmetry=False):
    if (not force_symmetry
        or crystal_symmetry is None
        or crystal_symmetry.space_group_info() is None):
      space_group_info = self.space_group_info()
    else:
      space_group_info = None
    return miller.set(
      crystal_symmetry=self.crystal_symmetry().join_symmetry(
        other_symmetry=crystal_symmetry,
        force=force_symmetry),
      indices=self.original_indices,
      anomalous_flag=True)

  def as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    if (base_array_info is None):
      base_array_info = miller.array_info(
        source_type="scalepack_no_merge_original_index")
    crystal_symmetry_from_file = self.crystal_symmetry()
    crystal_symmetry = crystal_symmetry_from_file.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=force_symmetry)
    result = miller.array(
      miller_set=self.unmerged_miller_set(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=True),
      data=self.i_obs,
      sigmas=self.sigmas)
    if (merge_equivalents):
      result = result.merge_equivalents().array()
    return (result
      .set_info(base_array_info.customized_copy(
        labels=["I(+)", "SIGI(+)", "I(-)", "SIGI(-)"],
        merged=merge_equivalents,
        crystal_symmetry_from_file=crystal_symmetry_from_file))
      .set_observation_type_xray_intensity())

  def as_miller_arrays(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None):
    return [self.as_miller_array(
      crystal_symmetry=crystal_symmetry,
      force_symmetry=force_symmetry,
      merge_equivalents=merge_equivalents,
      base_array_info=base_array_info)]

def combine_symops_and_symbol(space_group_from_ops, space_group_symbol):
  space_group_symbol = space_group_symbol.replace(" ","").upper()
  z = space_group_symbol[:1]
  if ("PABCIFRH".find(z) < 0):
    raise RuntimeError(
      "Cannot determine lattice centring type given space group symbol"
      " %s" % show_string(space_group_symbol))
  if (z == "P"):
    return sgtbx.space_group_info(group=space_group_from_ops)
  if (z == "H"):
    space_group_symbol = "R" + space_group_symbol[1:] + ":H"
    z = "R"
  elif (z == "R" and not space_group_symbol.endswith(":H")):
    if (space_group_symbol.endswith(":R")):
      z = None
    else:
      for s in space_group_from_ops:
        r_info = s.r().info()
        if (abs(r_info.type()) == 3):
          if (r_info.ev() == (0,0,1)):
            space_group_symbol = "R" + space_group_symbol[1:] + ":H"
            break
          elif (r_info.ev() == (1,1,1)):
            space_group_symbol += ":R"
            z = None
            break
  space_group_exp = sgtbx.space_group(space_group_from_ops)
  if (z is not None):
    try:
      space_group_exp.expand_conventional_centring_type(z)
    except RuntimeError:
      space_group_exp = None
  if (space_group_exp is not None):
    try:
      space_group_from_symbol = sgtbx.space_group_info(
        symbol=space_group_symbol).group()
    except RuntimeError:
      space_group_from_symbol = None
  if (   space_group_exp is None
      or space_group_from_symbol is None
      or space_group_exp != space_group_from_symbol):
    raise RuntimeError(
      "Symmetry operations in unmerged SCALEPACK file incompatible with"
      " space group symbol %s"
        % show_string(space_group_symbol))
  return sgtbx.space_group_info(group=space_group_exp)

def exercise_combine_symops_and_symbol():
  for symbols in sgtbx.space_group_symbol_iterator():
    space_group = sgtbx.space_group_info(
      symbol="Hall: %s" % symbols.hall()).group()
    space_group_from_ops = sgtbx.space_group()
    for s in list(space_group)[:space_group.order_p()]:
      space_group_from_ops.expand_smx(s)
    combined = combine_symops_and_symbol(
      space_group_from_ops=space_group_from_ops,
      space_group_symbol=symbols.universal_hermann_mauguin())
    assert combined.group() == space_group
    if (symbols.extension() in ["R", "H"]):
      combined = combine_symops_and_symbol(
        space_group_from_ops=space_group_from_ops,
        space_group_symbol=symbols.hermann_mauguin())
      assert combined.group() == space_group
      if (symbols.extension() == "H"):
        combined = combine_symops_and_symbol(
          space_group_from_ops=space_group_from_ops,
          space_group_symbol="H "+symbols.hermann_mauguin()[2:])
        assert combined.group() == space_group

def quick_test(file_name):
  from libtbx.utils import user_plus_sys_time
  t = user_plus_sys_time()
  s = reader(file_name)
  print "Time read:", t.delta()
  s.show_summary()
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
  m = s.as_miller_array(merge_equivalents=False).merge_equivalents()
  print "min redundancies:", flex.min(m.redundancies().data())
  print "max redundancies:", flex.max(m.redundancies().data())
  print "mean redundancies:", flex.mean(m.redundancies().data().as_double())
  s.as_miller_arrays()[0].show_summary()
  print

def run(args):
  exercise_combine_symops_and_symbol()
  file_names = """
p9/data/infl.sca
p9/data/peak.sca
p9/data/high.sca
rh-dehalogenase/data/auki_rd_1.sca
rh-dehalogenase/data/hgi2_rd_1.sca
rh-dehalogenase/data/hgki_rd_1.sca
rh-dehalogenase/data/ndac_rd_1.sca
rh-dehalogenase/data/rt_rd_1.sca
rh-dehalogenase/data/smac_1.sca
vmp/data/infl.sca
vmp/data/peak.sca
vmp/data/high.sca
bnl_2003/karen/shelxd/p123-unmerged.sca
bnl_2003/karen/shelxd/pk1-unmerged.sca
bnl_2003/karen/shelxd/pk12-unmerged.sca
bnl_2003/karen/shelxd/pk1234-unmerged.sca""".split()
  for file_name in file_names:
    for root_dir in args:
      fn = root_dir + "/" + file_name
      if (os.path.isfile(fn)):
        print "File name:", fn
        quick_test(fn)
  print "OK"

if (__name__ == "__main__"):
  run(sys.argv[1:])
