
from __future__ import division
from mmtbx.command_line import molprobity
from libtbx.easy_pickle import loads, dumps, dump
from libtbx.test_utils import show_diff
from libtbx.utils import null_out
import libtbx.load_env
from cStringIO import StringIO
import os

def exercise_molprobity () :
  if (not libtbx.env.has_module(name="probe")):
    print "Skipping exercise_molprobity(): probe not configured"
    return

  pdb_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3ifk.pdb",
    test=os.path.isfile)
  hkl_file = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/3ifk.mtz",
    test=os.path.isfile)
  if (pdb_file is None) :
    print "phenix_regression not available, skipping."
    return
  args1 = [
    pdb_file,
    "outliers_only=True",
  ]
  result = molprobity.run(args=args1, out=null_out())
  out1 = StringIO()
  result.show(out=out1)
  result = loads(dumps(result))
  out2 = StringIO()
  result.show(out=out2)
  assert (out2.getvalue() == out1.getvalue())
  dump("tst_molprobity.pkl", result)
  mc = result.as_multi_criterion_view()
  # percentiles
  out4 = StringIO()
  result.show_summary(out=out4, show_percentiles=True)
  assert ("""  Clashscore            =  49.96 (percentile: 1.0)""" in
    out4.getvalue())
  #result.show()
  assert (str(mc.data()[2]) == ' A   5  THR  rota,cb,clash')
  import mmtbx.validation.molprobity
  from iotbx import file_reader
  pdb_in = file_reader.any_file(pdb_file)
  hierarchy = pdb_in.file_object.construct_hierarchy()
  flags = mmtbx.validation.molprobity.molprobity_flags()
  flags.clashscore = False
  flags.model_stats = False
  flags.cbetadev = False
  result = mmtbx.validation.molprobity.molprobity(
    pdb_hierarchy=hierarchy,
    flags=flags)
  out3 = StringIO()
  result.show_summary(out=out3)
  assert not show_diff(out3.getvalue(), """\
  Ramachandran outliers =   1.76 %
                favored =  96.47 %
  Rotamer outliers      =  18.67 %
""")

if (__name__ == "__main__") :
  exercise_molprobity()
  print "OK"
