
from __future__ import absolute_import, division, print_function
from mmtbx.command_line import molprobity
import iotbx.pdb.hierarchy
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from cStringIO import StringIO
from six.moves import zip

# test for corner cases (synthetic data okay)
def exercise_synthetic():
  from mmtbx.regression import tst_build_alt_confs
  pdb_in = iotbx.pdb.hierarchy.input(pdb_string=tst_build_alt_confs.pdb_raw)
  xrs = pdb_in.input.xray_structure_simple()
  fc = abs(xrs.structure_factors(d_min=1.5).f_calc())
  flags = fc.resolution_filter(d_min=1.6).generate_r_free_flags()
  ls = fc.lone_set(other=flags)
  # case 1: no work set in high-res shell
  flags2 = ls.array(data=flex.bool(ls.size(), True))
  flags_all = flags.concatenate(other=flags2)
  mtz_out = fc.as_mtz_dataset(column_root_label="F")
  mtz_out.add_miller_array(flags_all, column_root_label="FreeR_flag")
  mtz_out.mtz_object().write("tst_molprobity_1.mtz")
  open("tst_molprobity_1.pdb", "w").write(tst_build_alt_confs.pdb_raw)
  args = [
    "tst_molprobity_1.pdb",
    "tst_molprobity_1.mtz",
    "--kinemage",
    "--maps",
    "flags.clashscore=False",
    "flags.xtriage=True",
  ]
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(out=out)
  # case 2: no test set in high-res shell
  flags2 = ls.array(data=flex.bool(ls.size(), False))
  flags_all = flags.concatenate(other=flags2)
  mtz_out = fc.as_mtz_dataset(column_root_label="F")
  mtz_out.add_miller_array(flags_all, column_root_label="FreeR_flag")
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(out=out)
  # case 3: multi-MODEL structure
  # XXX This is not a very sophisticated test - it only ensures that the
  # program does not crash.  We need a test for expected output...
  hierarchy = pdb_in.hierarchy
  model2 = hierarchy.only_model().detached_copy()
  hierarchy.append_model(model2)
  hierarchy.models()[0].id = "1"
  hierarchy.models()[1].id = "2"
  open("tst_molprobity_multi_model.pdb", "w").write(hierarchy.as_pdb_string())
  args = [
    "tst_molprobity_multi_model.pdb",
    "tst_molprobity_1.mtz",
    "--kinemage",
    "--maps",
  ]
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(out=out)
  # test rotamer distributions
  open("tst_molprobity_misc1.pdb", "w").write(tst_build_alt_confs.pdb_raw)
  args = [
    "tst_molprobity_1.pdb",
    "rotamer_library=8000",
  ]
  out = StringIO()
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  result.show(outliers_only=False, out=out)

def exercise_cdl():
  pdb_raw = """
ATOM   1270  N   LEU A 199       6.903  55.119  -0.416  1.00 25.48           N
ATOM   1271  CA  LEU A 199       7.726  56.192  -0.941  1.00 25.93           C
ATOM   1272  C   LEU A 199       6.996  56.972  -2.047  1.00 26.39           C
ATOM   1273  O   LEU A 199       7.020  58.180  -2.064  1.00 25.38           O
ATOM   1274  CB  LEU A 199       9.033  55.633  -1.490  1.00 25.66           C
ATOM   1278  N   ARG A 200       6.361  56.258  -2.980  1.00 27.31           N
ATOM   1279  CA  ARG A 200       5.576  56.913  -3.993  1.00 28.53           C
ATOM   1280  C   ARG A 200       4.520  57.823  -3.397  1.00 27.54           C
ATOM   1281  O   ARG A 200       4.397  58.949  -3.851  1.00 27.50           O
ATOM   1282  CB  ARG A 200       4.933  55.879  -4.899  1.00 30.38           C
ATOM   1289  N   ALA A 201       3.790  57.365  -2.357  1.00 26.90           N
ATOM   1290  CA  ALA A 201       2.764  58.200  -1.713  1.00 26.49           C
ATOM   1291  C   ALA A 201       3.406  59.407  -1.045  1.00 26.59           C
ATOM   1292  O   ALA A 201       2.866  60.516  -1.082  1.00 26.58           O
ATOM   1293  CB  ALA A 201       1.959  57.412  -0.715  1.00 25.11           C
ATOM   1294  N   ARG A 202       4.566  59.205  -0.419  1.00 25.66           N
ATOM   1295  CA  ARG A 202       5.245  60.296   0.240  1.00 26.93           C
ATOM   1296  C   ARG A 202       5.676  61.346  -0.767  1.00 26.93           C
ATOM   1297  O   ARG A 202       5.555  62.541  -0.489  1.00 25.79           O
ATOM   1298  CB  ARG A 202       6.493  59.779   0.996  1.00 28.25           C
ATOM   1305  N   ILE A 203       6.154  60.912  -1.931  1.00 26.99           N
ATOM   1306  CA  ILE A 203       6.611  61.848  -2.965  1.00 27.49           C
ATOM   1307  C   ILE A 203       5.430  62.674  -3.480  1.00 28.29           C
ATOM   1308  O   ILE A 203       5.548  63.905  -3.624  1.00 27.82           O
ATOM   1309  CB  ILE A 203       7.322  61.125  -4.075  1.00 28.09           C
ATOM   1313  N   SER A 204       4.288  62.025  -3.678  1.00 27.96           N
ATOM   1314  CA  SER A 204       3.119  62.736  -4.184  1.00 28.04           C
ATOM   1315  C   SER A 204       2.683  63.793  -3.199  1.00 28.16           C
ATOM   1316  O   SER A 204       2.311  64.910  -3.605  1.00 28.25           O
ATOM   1317  CB  SER A 204       1.962  61.780  -4.504  1.00 27.64           C
"""
  pdb_raw_2 = """\
REMARK   3    GEOSTD + MON.LIB. + CDL v1.2
""" + pdb_raw
  open("tst_molprobity_cdl_1.pdb", "w").write(pdb_raw)
  open("tst_molprobity_cdl_2.pdb", "w").write(pdb_raw_2)
  files = ["tst_molprobity_cdl_1.pdb","tst_molprobity_cdl_2.pdb"]
  rmsds = [0.9019, 0.8769]
  for file_name, rmsd, cdl_expected in zip(files, rmsds, [False, True]):
    result = molprobity.run(args=[file_name, "flags.clashscore=False"],
      ignore_missing_modules=True,
      out=null_out()).validation
    assert approx_equal(result.rms_angles(), rmsd, eps=0.001), rmsd
    if cdl_expected :
      out = StringIO()
      result.show(out=out)
      assert ("conformation-dependent library" in out.getvalue())
    else:
      out = StringIO()
      result.show(out=out)
      assert ("conformation-dependent library" not in out.getvalue())

if (__name__ == "__main__"):
  exercise_cdl()
  exercise_synthetic()
  print("OK")
