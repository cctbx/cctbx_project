
from __future__ import division
from __future__ import print_function
from cStringIO import StringIO
from libtbx.test_utils import show_diff
from libtbx.utils import null_out
from libtbx import easy_pickle
from libtbx import group_args

def exercise_heavy():
  from mmtbx.regression import make_fake_anomalous_data
  from mmtbx.command_line import validate_waters
  import mmtbx.ions.utils
  from iotbx.file_reader import any_file
  file_base = "tst_validate_waters_1"
  pdb_file = make_fake_anomalous_data.write_pdb_input_cd_cl(file_base=file_base)
  mtz_file = make_fake_anomalous_data.generate_mtz_file(
    file_base="tst_validate_waters_1",
    d_min=1.5,
    anomalous_scatterers=[
      group_args(selection="element CD", fp=-0.29, fdp=2.676),
      group_args(selection="element CL", fp=0.256, fdp=0.5),
    ])
  pdb_in = any_file(pdb_file)
  hierarchy = pdb_in.file_object.hierarchy
  hierarchy, n = mmtbx.ions.utils.anonymize_ions(hierarchy, log=null_out())
  hierarchy.write_pdb_file("%s_start.pdb" % file_base,
    crystal_symmetry=pdb_in.file_object.crystal_symmetry())
  args = ["tst_validate_waters_1_start.pdb", "tst_validate_waters_1.mtz",
    "skip_twin_detection=True"]
  results = validate_waters.run(args=args, out=null_out())
  out = StringIO()
  results.show(out=out)
  s = easy_pickle.dumps(results)
  r2 = easy_pickle.loads(s)
  out2 = StringIO()
  r2.show(out=out2)
  assert not show_diff(out.getvalue(), out2.getvalue())
  assert (results.n_bad >= 1) and (results.n_heavy == 2)
  # XXX statistics are approximate (probably a precision issue), so we can't
  # directly compare outputs

if (__name__ == "__main__"):
  exercise_heavy()
  print("OK")
