
from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from libtbx import easy_pickle
import libtbx.load_env
from six.moves import cStringIO as StringIO
import os
from six.moves import range

def exercise():
  for module in ["reduce", "probe", "phenix_regression"] :
    if (not libtbx.env.has_module(module)):
      print("%s not available, skipping" % module)
      return
  from mmtbx.command_line import validation_summary
  import iotbx.pdb.hierarchy
  import iotbx.pdb
  regression_pdb = libtbx.env.find_in_repositories(
    relative_path="phenix_regression/pdb/pdb1jxt.ent",
    test=os.path.isfile)
  out = StringIO()
  summary = validation_summary.run(args=[regression_pdb], out=out)
  assert approx_equal(summary.clashscore, 13.75, eps=0.001), \
    "clashscore %s is not 13.75(0.0001)" % summary.clashscore
  ss = easy_pickle.dumps(summary)
  sss = easy_pickle.loads(ss)
  out_1 = StringIO()
  out_2 = StringIO()
  summary.show(out=out_1)
  sss.show(out=out_2)
  assert out_1.getvalue() == out_2.getvalue()
  pdb_in = iotbx.pdb.input(regression_pdb)
  hierarchy = pdb_in.construct_hierarchy()
  new_hierarchy = iotbx.pdb.hierarchy.root()
  for i in range(5):
    model = hierarchy.only_model().detached_copy()
    model.id = str(i+1)
    new_hierarchy.append_model(model)
  import mmtbx.model
  model_object = mmtbx.model.manager(
          model_input=None,
          pdb_hierarchy = new_hierarchy,
          crystal_symmetry = pdb_in.crystal_symmetry())
  open("tst_validation_summary.pdb", "w").write(model_object.model_as_pdb())
  out2 = StringIO()
  summary = validation_summary.run(args=["tst_validation_summary.pdb"],
    out=out2)
  assert (type(summary).__name__ == 'ensemble')
  print("OK")

if (__name__ == "__main__"):
  exercise()
