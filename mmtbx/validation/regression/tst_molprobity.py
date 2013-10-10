
from __future__ import division
from mmtbx.command_line import molprobity
from libtbx.easy_pickle import loads, dumps, dump
from libtbx.utils import null_out
import libtbx.load_env
from cStringIO import StringIO
import os

def exercise_molprobity () :
  pdb_in = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/3ifk.pdb",
    test=os.path.isfile)
  hkl_in = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/reflection_files/3ifk.mtz",
    test=os.path.isfile)
  if (pdb_in is None) :
    print "phenix_regression not available, skipping."
    return
  args1 = [
    pdb_in,
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
  result.show()
  assert (str(mc.data()[2]) == ' A   5  THR  rota,cb,clash')
  print "OK"

if (__name__ == "__main__") :
  exercise_molprobity()
