
from __future__ import division
from mmtbx.command_line import molprobity
import mmtbx.validation.molprobity
import iotbx.pdb.hierarchy
from scitbx.array_family import flex
from libtbx.utils import null_out
from cStringIO import StringIO

# test for corner cases (synthetic data okay)
def exercise_synthetic () :
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
    "flags.clashscore=False",
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
    "rotamer_library=500",
    "outliers_only=False",
    "flags.clashscore=False",
  ]
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  out = StringIO()
  result.show(outliers_only=False, out=out)
  assert ("""A   7  TYR              m-85  91.83   299.6,92.0""" in
          out.getvalue())
  args2 = [
    "tst_molprobity_1.pdb",
    "rotamer_library=8000",
    "outliers_only=False",
    "flags.clashscore=False",
  ]
  out = StringIO()
  result = molprobity.run(args=args2,
    ignore_missing_modules=True,
    out=null_out()).validation
  result.show(outliers_only=False, out=out)
  assert ("""   A   7  TYR              m-85  93.10   299.6,92.0""" in
    out.getvalue()), out.getvalue()

if (__name__ == "__main__") :
  exercise_synthetic()
  print "OK"
