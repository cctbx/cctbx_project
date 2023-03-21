
# TODO more tests

from __future__ import absolute_import, division, print_function
from iotbx.pdb import hierarchy
from mmtbx.command_line import molprobity
import iotbx.pdb
from scitbx.array_family import flex
from libtbx.utils import null_out
import libtbx.load_env
import random

# Deuterium as ligand - formerly crashed real-space correlation
def exercise_01():
  pdb_raw = """\
CRYST1   10.000   15.000   10.000  90.00  90.00  90.00 P 1
ATOM   6407  N   GLY A 388      -0.783   9.368 -16.436  1.00 51.96           N
ATOM   6408  CA  GLY A 388      -0.227   9.888 -15.197  1.00 54.04           C
ATOM   6409  C   GLY A 388      -0.637  11.320 -14.897  1.00 55.86           C
ATOM   6410  O   GLY A 388      -1.728  11.738 -15.347  1.00 56.70           O
ATOM   6411  OXT GLY A 388       0.129  12.024 -14.203  1.00 56.98           O
ATOM   6412  D   GLY A 388      -0.460   9.727 -17.309  1.00 51.44           D
ATOM   6413  HA2 GLY A 388      -0.561   9.258 -14.385  1.00 54.07           H
ATOM   6414  HA3 GLY A 388       0.843   9.835 -15.243  1.00 54.13           H
TER    6415      GLY A 388
HETATM 6416  D   D8U A 401       2.236   5.695 -12.992  1.00 15.23           D
HETATM 6417  O   DOD A1001      -4.151   4.107 -16.592  1.00 13.40           O
HETATM 6418  D1  DOD A1001      -4.760   3.026 -11.326  1.00 15.45           D
HETATM 6419  D2  DOD A1001      -4.625   2.741 -13.845  1.00 14.81           D
"""
  random.seed(12345)
  flex.set_random_seed(12345)
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_raw)
  xrs = pdb_in.xray_structure_simple()
  hierarchy = pdb_in.construct_hierarchy()
  hierarchy.atoms().reset_i_seq()
  open("tst_validate_experimental.pdb", "w").write(
    hierarchy.as_pdb_string(crystal_symmetry=xrs))
  f_calc =abs( xrs.structure_factors(d_min=2.0).f_calc())
  f_calc.set_observation_type_xray_amplitude()
  flags = f_calc.generate_r_free_flags()
  mtz = f_calc.as_mtz_dataset(column_root_label="F")
  mtz.add_miller_array(flags, column_root_label="FreeR_flag")
  mtz.mtz_object().write("tst_validate_experimental.mtz")
  args = [
    "tst_validate_experimental.pdb",
    "tst_validate_experimental.mtz",
  ]
  args.append("flags.clashscore=%s" % libtbx.env.has_module("probe"))
  result = molprobity.run(args=args,
    ignore_missing_modules=True,
    out=null_out()).validation
  assert result.real_space is not None
  ds = result.data_stats
  assert (ds.n_free > 0)

if (__name__ == "__main__"):
  exercise_01()
  print("OK")
