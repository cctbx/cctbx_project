from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
import mmtbx.model
from libtbx.utils import null_out
import iotbx.pdb
from mmtbx.validation.rama_z import rama_z
import os

fname = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/secondary_structure/regularize_from_pdb_lib/lib8.pdb",
    test=os.path.isfile)

def check_function():
  inp = iotbx.pdb.input(fname)
  model = mmtbx.model.manager(model_input=inp)
  zs = rama_z(model, log=null_out())
  z_scores = zs.get_z_scores()
  ss_cont = zs.get_residue_counts()
  # print (z_scores)
  # print (ss_cont)
  expected_z =  {'H': -2.7111077448509953, 'S': None, 'L': -0.7253113914835396, 'W': -1.0306993840276084}
  expeted_ss = {'H': 277, 'S': 0, 'L': 15041, 'W': 15318}
  for k in expected_z:
    if z_scores[k] is not None:
      assert approx_equal( z_scores[k][0], expected_z[k])
      if k == 'W':
        assert 0.01 < z_scores[k][1] < 0.1, z_scores[k][1]
    if k != 'weighted_mean':
      assert approx_equal( ss_cont[k], expeted_ss[k] )

def check_cmd_line():
  cmd = "mmtbx.rama_z %s" % fname
  r = easy_run.fully_buffered(cmd)
  stdout = r.stdout_lines
  # print ("\n".join(stdout))
  expected_strs = [
      ["z-score whole: -1.031", "residues: 15318"],
      ["z-score helix: -2.711", "residues: 277"],
      ["z-score sheet: None,", "residues: 0"],
      ["z-score loop : -0.725", "residues: 15041"]]
  for res, expected in zip(stdout[-9:-6], expected_strs):
    for exp_l in expected:
      assert res.find(exp_l) >= 0, "res: '%s', exp: '%s'" % (res, exp_l)

if __name__ == '__main__':
  check_function()
  check_cmd_line()
  print("OK")
