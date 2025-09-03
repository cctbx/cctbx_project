from __future__ import absolute_import, division, print_function
from libtbx import easy_run
from libtbx.test_utils import approx_equal, assert_lines_in_text
import mmtbx.model
from libtbx.utils import null_out
import iotbx.pdb
from mmtbx.validation.rama_z import rama_z

import mmtbx
from pathlib import Path
data_dir = Path(mmtbx.__file__).parent / 'regression' / 'pdbs'
fname = str( data_dir / 'p9.pdb')


cryst1 = "CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1          \n"

def check_function():
  inp = iotbx.pdb.input(fname)
  model = mmtbx.model.manager(model_input=inp)
  zs = rama_z([model], log=null_out())
  z_scores = zs.get_z_scores()
  ss_cont = zs.get_residue_counts()
  # print (z_scores)
  # print (ss_cont)
  expected_z =  {'H': None, 'S': (-0.057428666470734, 0.5840017477579902),
      'L': (-0.3588028726184504, 0.6941226745661744),
      'W': (-0.4019606027769244, 0.6621289642029733)}
  expeted_ss = {'H': 0, 'S': 63, 'L': 71, 'W': 134}
  for k in expected_z:
    if z_scores[k] is not None:
      assert approx_equal( z_scores[k], expected_z[k], eps=1e-5)
      assert approx_equal( ss_cont[k], expeted_ss[k] )
  # check how separate scores translate to whole
  s_score = (z_scores['S'][0] * zs.calibration_values['S'][1] + zs.calibration_values['S'][0]) * ss_cont['S']
  l_score = (z_scores['L'][0] * zs.calibration_values['L'][1] + zs.calibration_values['L'][0]) * ss_cont['L']
  w_score = ((s_score + l_score)/(ss_cont['S']+ss_cont['L']) - zs.calibration_values['W'][0]) / zs.calibration_values['W'][1]
  # print ("reconstructed:", w_score, z_scores['W'][0])
  assert approx_equal(w_score, z_scores['W'][0])

def check_cmd_line():
  cmd = "mmtbx.rama_z %s" % fname
  r = easy_run.fully_buffered(cmd)
  assert r.return_code == 0
  stdout = r.stdout_lines
  # print ("\n".join(stdout))
  assert_lines_in_text("\n".join(stdout), """\
      whole: -0.40 (0.66), residues: 134
      helix:  None (None), residues: 0
      sheet: -0.06 (0.58), residues: 63
      loop : -0.36 (0.69), residues: 71""")

def check_cmd_line_cryst1(prefix="tst_rama_z_01_cryst1"):
  with open(fname, 'r') as f:
    pdbtext = f.read()
  with open(prefix+'.pdb', 'w') as f:
    f.write(cryst1)
    f.write(pdbtext)
  cmd = "mmtbx.rama_z %s" % (prefix+'.pdb')
  r = easy_run.fully_buffered(cmd)
  assert r.return_code == 0
  stdout = r.stdout_lines
  # print ("\n".join(stdout))
  assert_lines_in_text("\n".join(stdout), """\
      whole: -0.40 (0.66), residues: 134
      helix:  None (None), residues: 0
      sheet: -0.06 (0.58), residues: 63
      loop : -0.36 (0.69), residues: 71""")


if __name__ == '__main__':
  check_function()
  check_cmd_line()
  check_cmd_line_cryst1()
  print("OK")
