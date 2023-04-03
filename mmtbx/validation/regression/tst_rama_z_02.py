from __future__ import absolute_import, division, print_function
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import assert_lines_in_text
import os

fname = libtbx.env.find_in_repositories(
    relative_path="cctbx_project/mmtbx/regression/pdbs/1ucs_cutted_xyz_rounded.pdb",
    test=os.path.isfile)

def check_cmd_line():
  cmd = "mmtbx.rama_z %s" % fname
  r = easy_run.fully_buffered(cmd)
  stdout = r.stdout_lines
  # stderr = r.stderr_lines
  # print ("\n".join(stdout))
  # print ("\n".join(stderr))
  # assert r.return_code == 0
  assert_lines_in_text("\n".join(stdout), """\
      whole: -7.73 (0.57), residues: 38
      helix: -5.87 (0.32), residues: 6
      sheet:  None (None), residues: 0
      loop : -5.60 (0.49), residues: 32""")

if __name__ == '__main__':
  check_cmd_line()
  print("OK")
