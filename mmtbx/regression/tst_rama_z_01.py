from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal, assert_lines_in_text
import mmtbx.model
from libtbx.utils import null_out
import iotbx.pdb
from mmtbx.validation.rama_z import rama_z
import os

fname = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/regression/pdbs/p9.pdb",
    test=os.path.isfile)

def check_function():
  inp = iotbx.pdb.input(fname)
  model = mmtbx.model.manager(model_input=inp)
  zs = rama_z(model, params=None, log=null_out())
  z_scores = zs.get_z_scores()
  ss_cont = zs.get_residue_counts()
  # print (z_scores)
  # print (ss_cont)
  expected_z =  {'H': None, 'S': -0.057428666470734, 'L': -0.3588028726184504,
      'W': -0.4019606027769244}
  expeted_ss = {'H': 0, 'S': 63, 'L': 71, 'W': 134}
  for k in expected_z:
    if z_scores[k] is not None:
      assert approx_equal( z_scores[k], expected_z[k] )
    if k != 'weighted_mean':
      assert approx_equal( ss_cont[k], expeted_ss[k] )

def check_cmd_line():
  cmd = "mmtbx.rama_z %s" % fname
  r = easy_run.fully_buffered(cmd)
  stdout = r.stdout_lines
  # print ("\n".join(stdout))

  assert_lines_in_text("\n".join(stdout), """\
      z-score whole: -0.402, residues: 134
      z-score helix: None, residues: 0
      z-score sheet: -0.057, residues: 63
      z-score loop : -0.359, residues: 71
      """)

if __name__ == '__main__':
  check_function()
  check_cmd_line()
  print("OK")
