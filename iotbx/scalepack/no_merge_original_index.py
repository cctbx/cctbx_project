from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
from cctbx.array_family import flex
from libtbx.str_utils import show_string
import warnings
import sys, os

import boost_adaptbx.boost.python as bp
from six.moves import range
scalepack_ext = bp.import_ext("iotbx_scalepack_ext")

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

def writer(
    i_obs,
    file_object=None,
    file_name=None,
    batch_numbers=None,
    spindle_flags=None,
    scale_intensities_for_scalepack_merge=False,
    out=sys.stdout):
  n_refl = len(i_obs.indices())
  #assert (not i_obs.is_unique_set_under_symmetry()) # TT 2014-01-12 not necessary?
  assert i_obs.is_xray_intensity_array() and i_obs.sigmas() is not None
  # create substitute for batch number and spindle flag - these are impossible
  # to determine from the input array
  if (type(batch_numbers).__name__ == "array"):
    assert batch_numbers.indices().all_eq(i_obs.indices())
    batch_numbers = batch_numbers.data()
  elif (batch_numbers is None):
    batch_numbers = flex.int(n_refl, 0)
  if (type(spindle_flags).__name__ == "array"):
    assert spindle_flags.indices().all_eq(i_obs.indices())
    spindle_flags = spindle_flags.data()
  elif (spindle_flags is None):
    spindle_flags = flex.int(n_refl, 0)
  assert len(batch_numbers) == len(spindle_flags) == n_refl
  assert ([file_object, file_name].count(None) == 1)
  if (file_object is None):
    file_object = open(file_name, "w")
  space_group = i_obs.space_group()
  # generate array for the centric/I+/I- flag
  centric_flags = i_obs.centric_flags().data()
  i_obs_asu = i_obs.as_non_anomalous_array().map_to_asu()
  i_obs_asu_anom = i_obs.map_to_asu()
  friedel_mate_flags = ~(i_obs_asu.indices() == i_obs_asu_anom.indices())
  centric_tags = flex.int(n_refl, 1)
  centric_tags.set_selected(friedel_mate_flags.iselection(), 2)
  centric_tags.set_selected(centric_flags.iselection(), 0)
  # generate isym array
  uniq_indices = i_obs.indices().deep_copy()
  #uniq_indices_anom = i_obs.indices().deep_copy()
  isym = flex.int(n_refl, 0)
  miller.map_to_asu_isym(space_group.type(), False, uniq_indices, isym)
  #miller.map_to_asu_isym(space_group.type(), False, uniq_indices_anom,
  #  isym)
  # write out symmetry operators
  n_smx = space_group.n_smx()
  file_object.write("%5d %s\n" % (n_smx, i_obs.space_group_info()))
  for smx in space_group.smx():
    smx_rot = smx.as_int_array()[0:9]
    smx_tra = smx.as_int_array()[9:]
    for r in smx_rot :
      file_object.write("%3d" % r)
    file_object.write("\n")
    for t in smx_tra :
      file_object.write("%3d" % t)
    file_object.write("\n")
  # write out reflections
  if scale_intensities_for_scalepack_merge: # 2014-01-07 TT
    from iotbx.scalepack.merge import scale_intensities_if_necessary
    i_obs=scale_intensities_if_necessary(i_obs,out=out)

  from iotbx.scalepack.merge import format_f8_1_or_i8 # Sorry if out of range
  for i_refl, (h,k,l) in enumerate(i_obs.indices()):
    (h_asu, k_asu, l_asu) = uniq_indices[i_refl]
    c_flag = centric_tags[i_refl]
    asu_number = abs(isym[i_refl]) + 1 # XXX is this correct?
    spindle_flag = spindle_flags[i_refl]
    batch = batch_numbers[i_refl]
    i_formatted=format_f8_1_or_i8((h,k,l),"intensity",i_obs.data()[i_refl])
    s_formatted=format_f8_1_or_i8((h,k,l),"sigma",i_obs.sigmas()[i_refl])

    file_object.write("%4d%4d%4d%4d%4d%4d%6d%2d%2d%3d%s%s\n" %
      (h, k, l, h_asu, k_asu, l_asu, batch, c_flag, spindle_flag,
      asu_number, i_formatted,s_formatted,))
  file_object.close()

class reader(object):
  def __init__(self, file_name, header_only=False):
    self.file_name = os.path.normpath(file_name)
    with open(file_name) as f:
      line = f.readline()
      assert line[5] == " "
      n_sym_ops_from_file = int(line[:5].strip())
      assert n_sym_ops_from_file > 0
      self.space_group_symbol = line[6:].strip()
      self.space_group_from_ops = sgtbx.space_group()
      for i in range(n_sym_ops_from_file):
        line = f.readline().rstrip()
        assert len(line) == 27
        r = sgtbx.rot_mx([int(line[j*3:(j+1)*3]) for j in range(9)], 1)
        line = f.readline().rstrip()
        assert len(line) == 9
        t = sgtbx.tr_vec([int(line[j*3:(j+1)*3]) for j in range(3)], 12)
        self.space_group_from_ops.expand_smx(sgtbx.rt_mx(r, t))
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
    print(prefix + "File name:", show_string(self.file_name), file=out)
    print(prefix + "Space group symbol:", \
      show_string(self.space_group_symbol), file=out)
    try: space_group_info = self.space_group_info()
    except KeyboardInterrupt: raise
    except Exception: pass
    else:
      space_group_info.show_summary(
        f=out, prefix=prefix+"Space group from operations: ")
    if (self.original_indices is not None):
      print(prefix + "Number of original indices:", \
        self.original_indices.size(), file=out)

  def crystal_symmetry(self):
    return crystal.symmetry(
      unit_cell=None,
      space_group_info=self.space_group_info())

  def unmerged_miller_set(self,
                          crystal_symmetry=None,
                          force_symmetry=False,
                          anomalous=True):
    if anomalous is None:
      anomalous = True
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
      anomalous_flag=anomalous)

  def as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        merge_equivalents=True,
        base_array_info=None,
        anomalous=True):
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
        force_symmetry=True,
        anomalous=anomalous),
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
        base_array_info=None,
        anomalous=True):
    return [
      self.as_miller_array(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=force_symmetry,
        merge_equivalents=merge_equivalents,
        base_array_info=base_array_info,
        anomalous=anomalous,
      ),
      self.batch_as_miller_array(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=force_symmetry,
        base_array_info=base_array_info,
        anomalous=anomalous,
      ),
    ]

  def batch_as_miller_array(self,
        crystal_symmetry=None,
        force_symmetry=False,
        base_array_info=None,
        anomalous=True):
    if (base_array_info is None):
      base_array_info = miller.array_info(
        source_type="scalepack_no_merge_original_index")
    crystal_symmetry_from_file = self.crystal_symmetry()
    crystal_symmetry = crystal_symmetry_from_file.join_symmetry(
      other_symmetry=crystal_symmetry,
      force=force_symmetry)
    return miller.array(
      miller_set=self.unmerged_miller_set(
        crystal_symmetry=crystal_symmetry,
        force_symmetry=True,
        anomalous=anomalous),
      data=self.batch_numbers).set_info(
        base_array_info.customized_copy(
          labels=["BATCH"],
          crystal_symmetry_from_file=crystal_symmetry_from_file))

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
    if space_group_from_symbol:
      warnings.warn("""
WARNING:
  Symmetry operations in input file are for space group %(space_group_exp)s
  However space group symbol is: %(space_group_symbol)s
  This may be a format error in the Scalepack file!
  Using %(space_group_symbol)s
""" % {"space_group_exp" : str(space_group_exp.info()),
       "space_group_symbol" : show_string(space_group_symbol), },
        UserWarning, stacklevel=10)
      space_group_exp = space_group_from_symbol
    else:
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
  print("Time read:", t.delta())
  s.show_summary()
  print(tuple(s.original_indices[:3]))
  print(tuple(s.unique_indices[:3]))
  print(tuple(s.batch_numbers[:3]))
  print(tuple(s.centric_tags[:3]))
  print(tuple(s.spindle_flags[:3]))
  print(tuple(s.asymmetric_unit_indices[:3]))
  print(tuple(s.i_obs[:3]))
  print(tuple(s.sigmas[:3]))
  print(tuple(s.original_indices[-3:]))
  print(tuple(s.unique_indices[-3:]))
  print(tuple(s.batch_numbers[-3:]))
  print(tuple(s.centric_tags[-3:]))
  print(tuple(s.spindle_flags[-3:]))
  print(tuple(s.asymmetric_unit_indices[-3:]))
  print(tuple(s.i_obs[-3:]))
  print(tuple(s.sigmas[-3:]))
  m = s.as_miller_array(merge_equivalents=False).merge_equivalents()
  print("min redundancies:", flex.min(m.redundancies().data()))
  print("max redundancies:", flex.max(m.redundancies().data()))
  print("mean redundancies:", flex.mean(m.redundancies().data().as_double()))
  s.as_miller_arrays()[0].show_summary()
  print()

def run(args):
  exercise_combine_symops_and_symbol()

  file_names = """
p9/infl.sca
p9/peak.sca
p9/high.sca
rh-dehalogenase/auki_rd_1.sca
rh-dehalogenase/hgi2_rd_1.sca
rh-dehalogenase/hgki_rd_1.sca
rh-dehalogenase/ndac_rd_1.sca
rh-dehalogenase/rt_rd_1.sca
rh-dehalogenase/smac_1.sca
vmp/infl.sca
vmp/peak.sca
vmp/high.sca""".split()
  for file_name in file_names:
    for root_dir in args:
      fn = root_dir + "/" + file_name
      if (os.path.isfile(fn)):
        print("File name:", fn)
        quick_test(fn)
  print("OK")

if (__name__ == "__main__"):
  run(sys.argv[1:])
