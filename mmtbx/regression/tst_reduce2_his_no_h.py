"""reduce2 with flip Movers: a His keeps at least one ring H. MoverHisFlip
offered a state with both ring H removed (preference -1.0); in 4KW8 His181
picked it (NE2 near Thr62 OG1, ND1 near Ala179 O), leaving a neutral ring with
no H on either N. Model: 4KW8 Thr62, Ala179-His181 (Asp180 as deposited,
without OD1/OD2) and HOH 419, the water Thr62 HG1 can reach.

A His locked by set_flip_states (one orientation, 3 states) is recorded as
flipped when its Mover is added, never from its final state. Models: 6VW1
His B505 (locked Unflipped) and 3J4P His A428 (locked Flipped), both ending
with HD1 only, the state at index len//2 of a one-orientation Mover.
"""
from __future__ import absolute_import, division, print_function
import os
import iotbx.pdb
import mmtbx.model
from libtbx.utils import null_out
from iotbx.cli_parser import run_program
from mmtbx.programs import reduce2

pdb_str = """\
CRYST1   51.623   63.021   66.347  90.00  90.00  90.00 P 21 21 21    4
ATOM    837  N   THR A  62     102.257  -0.330  66.504  1.00 15.42           N
ATOM    838  CA  THR A  62     102.127   0.059  67.898  1.00 15.93           C
ATOM    839  C   THR A  62     103.493   0.431  68.504  1.00 15.92           C
ATOM    840  O   THR A  62     103.567   1.248  69.406  1.00 16.16           O
ATOM    841  CB  THR A  62     101.426  -1.023  68.755  1.00 17.18           C
ATOM    842  OG1 THR A  62     102.128  -2.272  68.653  1.00 17.70           O
ATOM    843  CG2 THR A  62      99.965  -1.194  68.324  1.00 17.55           C
ATOM   2538  N   ALA A 179      97.566 -10.157  68.988  1.00 24.17           N
ATOM   2539  CA  ALA A 179      98.749  -9.722  69.660  1.00 26.95           C
ATOM   2540  C   ALA A 179      98.335  -9.088  70.974  1.00 30.34           C
ATOM   2541  O   ALA A 179      97.515  -8.181  70.990  1.00 33.89           O
ATOM   2542  CB  ALA A 179      99.515  -8.763  68.791  1.00 27.13           C
ATOM   2548  N   ASP A 180      98.898  -9.605  72.067  1.00 31.16           N
ATOM   2549  CA  ASP A 180      98.674  -9.112  73.423  1.00 22.43           C
ATOM   2550  C   ASP A 180      99.846  -8.238  73.884  1.00 21.40           C
ATOM   2551  O   ASP A 180     100.952  -8.733  74.065  1.00 28.44           O
ATOM   2552  CB  ASP A 180      98.560 -10.302  74.394  1.00 24.21           C
ATOM   2553  CG  ASP A 180      97.276 -11.095  74.226  1.00 25.52           C
ATOM   2558  N   HIS A 181      99.592  -6.956  74.093  1.00 24.57           N
ATOM   2559  CA  HIS A 181     100.611  -5.985  74.508  1.00 23.83           C
ATOM   2560  C   HIS A 181     100.559  -5.740  76.007  1.00 26.15           C
ATOM   2561  O   HIS A 181      99.487  -5.585  76.578  1.00 23.51           O
ATOM   2562  CB  HIS A 181     100.379  -4.640  73.826  1.00 20.92           C
ATOM   2563  CG  HIS A 181     100.322  -4.723  72.338  1.00 24.35           C
ATOM   2564  ND1 HIS A 181      99.569  -5.664  71.669  1.00 26.53           N
ATOM   2565  CD2 HIS A 181     100.917  -3.969  71.386  1.00 24.37           C
ATOM   2566  CE1 HIS A 181      99.720  -5.499  70.369  1.00 24.79           C
ATOM   2567  NE2 HIS A 181     100.523  -4.468  70.170  1.00 25.17           N
HETATM 3296  O   HOH A 419     104.349  -3.083  67.365  1.00  5.42           O
END
"""

def test_his181_keeps_ring_h():
  prefix = "tst_reduce2_his_no_h"
  fn_in, fn_out, fn_desc = prefix + ".pdb", prefix + "_out.pdb", prefix + "_out.txt"
  with open(fn_in, "w") as f:
    f.write(pdb_str)
  for fn in (fn_out, fn_desc):
    if os.path.exists(fn): os.remove(fn)
  run_program(program_class=reduce2.Program, logger=null_out(), args=[fn_in,
    "output.filename=" + fn_out, "output.description_file_name=" + fn_desc,
    "overwrite=True", "add_flip_movers=True"])
  with open(fn_out) as f:
    model_out = f.read()
  with open(fn_desc) as f:
    desc = f.read()
  poses = [l for l in desc.splitlines() if "HisFlip at chain A HIS 181" in l]
  assert len(poses) == 1, desc
  assert "HD1NotPlaced HE2NotPlaced" not in poses[0], poses[0]
  ring_h = [l[12:16].strip() for l in model_out.splitlines()
    if l.startswith("ATOM") and l[17:26] == "HIS A 181" and l[12:16] in (" HD1", " HE2")]
  assert len(ring_h) >= 1, model_out

pdb_6vw1_his505 = """\
CRYST1   80.435  118.034  112.075  90.00  93.12  90.00 P 1 21 1
ATOM  13550  N   ARG B 273      56.889  27.909 122.225  1.00 -1.00           N
ATOM  13551  CA  ARG B 273      58.147  27.466 121.637  1.00 -1.00           C
ATOM  13552  C   ARG B 273      59.162  27.093 122.707  1.00 -1.00           C
ATOM  13553  O   ARG B 273      59.951  26.160 122.516  1.00 -1.00           O
ATOM  13554  CB  ARG B 273      58.711  28.558 120.726  1.00 -1.00           C
ATOM  13555  CG  ARG B 273      60.002  28.172 120.022  1.00 -1.00           C
ATOM  13556  CD  ARG B 273      60.565  29.327 119.211  1.00 -1.00           C
ATOM  13557  NE  ARG B 273      61.882  29.010 118.665  1.00 -1.00           N
ATOM  13558  CZ  ARG B 273      62.725  29.913 118.176  1.00 -1.00           C
ATOM  13559  NH1 ARG B 273      62.393  31.197 118.165  1.00 -1.00           N
ATOM  13560  NH2 ARG B 273      63.903  29.532 117.701  1.00 -1.00           N
ATOM  17229  N   PHE B 504      57.300  34.456 114.078  1.00 -1.00           N
ATOM  17230  CA  PHE B 504      58.198  35.260 113.256  1.00 -1.00           C
ATOM  17231  C   PHE B 504      58.781  34.454 112.104  1.00 -1.00           C
ATOM  17232  O   PHE B 504      58.761  34.900 110.953  1.00 -1.00           O
ATOM  17233  CB  PHE B 504      59.322  35.840 114.115  1.00 -1.00           C
ATOM  17234  CG  PHE B 504      60.460  36.404 113.313  1.00 -1.00           C
ATOM  17235  CD1 PHE B 504      60.373  37.667 112.752  1.00 -1.00           C
ATOM  17236  CD2 PHE B 504      61.614  35.665 113.112  1.00 -1.00           C
ATOM  17237  CE1 PHE B 504      61.418  38.182 112.007  1.00 -1.00           C
ATOM  17238  CE2 PHE B 504      62.660  36.173 112.368  1.00 -1.00           C
ATOM  17239  CZ  PHE B 504      62.562  37.433 111.816  1.00 -1.00           C
ATOM  17249  N   HIS B 505      59.338  33.277 112.400  1.00 -1.00           N
ATOM  17250  CA  HIS B 505      60.006  32.500 111.361  1.00 -1.00           C
ATOM  17251  C   HIS B 505      59.045  32.099 110.251  1.00 -1.00           C
ATOM  17252  O   HIS B 505      59.439  32.032 109.081  1.00 -1.00           O
ATOM  17253  CB  HIS B 505      60.659  31.259 111.970  1.00 -1.00           C
ATOM  17254  CG  HIS B 505      61.669  31.565 113.030  1.00 -1.00           C
ATOM  17255  ND1 HIS B 505      61.325  32.098 114.254  1.00 -1.00           N
ATOM  17256  CD2 HIS B 505      63.015  31.421 113.048  1.00 -1.00           C
ATOM  17257  CE1 HIS B 505      62.415  32.264 114.982  1.00 -1.00           C
ATOM  17258  NE2 HIS B 505      63.455  31.862 114.273  1.00 -1.00           N
ATOM  17352  N   PHE B 512      62.797  26.230 109.850  1.00 -1.00           N
ATOM  17353  CA  PHE B 512      62.710  26.002 111.285  1.00 -1.00           C
ATOM  17354  C   PHE B 512      63.381  24.725 111.775  1.00 -1.00           C
ATOM  17355  O   PHE B 512      63.738  24.662 112.956  1.00 -1.00           O
ATOM  17356  CB  PHE B 512      61.242  25.994 111.738  1.00 -1.00           C
ATOM  17357  CG  PHE B 512      61.040  26.506 113.132  1.00 -1.00           C
ATOM  17358  CD1 PHE B 512      61.287  27.834 113.437  1.00 -1.00           C
ATOM  17359  CD2 PHE B 512      60.590  25.666 114.135  1.00 -1.00           C
ATOM  17360  CE1 PHE B 512      61.102  28.312 114.720  1.00 -1.00           C
ATOM  17361  CE2 PHE B 512      60.402  26.136 115.421  1.00 -1.00           C
ATOM  17362  CZ  PHE B 512      60.658  27.461 115.714  1.00 -1.00           C
ATOM  17415  N   TYR B 515      66.574  24.376 114.205  1.00 -1.00           N
ATOM  17416  CA  TYR B 515      66.654  24.361 115.662  1.00 -1.00           C
ATOM  17417  C   TYR B 515      66.687  22.936 116.202  1.00 -1.00           C
ATOM  17418  O   TYR B 515      67.436  22.640 117.139  1.00 -1.00           O
ATOM  17419  CB  TYR B 515      65.501  25.172 116.248  1.00 -1.00           C
ATOM  17420  CG  TYR B 515      65.589  26.629 115.852  1.00 -1.00           C
ATOM  17421  CD1 TYR B 515      66.587  27.444 116.372  1.00 -1.00           C
ATOM  17422  CD2 TYR B 515      64.705  27.181 114.935  1.00 -1.00           C
ATOM  17423  CE1 TYR B 515      66.689  28.772 116.011  1.00 -1.00           C
ATOM  17424  CE2 TYR B 515      64.799  28.513 114.565  1.00 -1.00           C
ATOM  17425  CZ  TYR B 515      65.795  29.302 115.107  1.00 -1.00           C
ATOM  17426  OH  TYR B 515      65.899  30.628 114.748  1.00 -1.00           O
END
"""

pdb_3j4p_his428 = """\
ATOM   1670  N   TYR A 426       0.848  14.209 100.340  1.00 19.22           N
ATOM   1671  CA  TYR A 426       2.067  13.630 100.868  1.00 22.74           C
ATOM   1672  C   TYR A 426       2.650  14.376 102.054  1.00 27.50           C
ATOM   1673  O   TYR A 426       2.366  15.555 102.258  1.00 14.15           O
ATOM   1674  CB  TYR A 426       3.104  13.500  99.753  1.00 24.86           C
ATOM   1675  CG  TYR A 426       3.358  14.762  98.954  1.00 18.61           C
ATOM   1676  CD1 TYR A 426       4.087  15.825  99.493  1.00 25.45           C
ATOM   1677  CD2 TYR A 426       2.918  14.872  97.631  1.00 26.45           C
ATOM   1678  CE1 TYR A 426       4.382  16.965  98.729  1.00 18.86           C
ATOM   1679  CE2 TYR A 426       3.206  16.007  96.860  1.00 29.42           C
ATOM   1680  CZ  TYR A 426       3.940  17.047  97.415  1.00 27.08           C
ATOM   1681  OH  TYR A 426       4.244  18.156  96.655  1.00 15.84           O
ATOM   1687  N   HIS A 428       6.121  15.617 104.171  1.00 13.01           N
ATOM   1688  CA  HIS A 428       7.500  15.956 103.885  1.00 15.97           C
ATOM   1689  C   HIS A 428       8.519  15.210 104.722  1.00 24.32           C
ATOM   1690  O   HIS A 428       8.222  14.750 105.822  1.00 15.38           O
ATOM   1691  CB  HIS A 428       7.681  17.460 104.007  1.00 14.02           C
ATOM   1692  CG  HIS A 428       7.011  18.214 102.906  1.00 15.60           C
ATOM   1693  ND1 HIS A 428       5.683  18.028 102.587  1.00 27.83           N
ATOM   1694  CD2 HIS A 428       7.495  19.101 102.005  1.00  9.02           C
ATOM   1695  CE1 HIS A 428       5.380  18.764 101.533  1.00 26.05           C
ATOM   1696  NE2 HIS A 428       6.462  19.425 101.160  1.00 21.50           N
ATOM   3140  N   GLN A 609       1.256  24.384 102.485  1.00  8.17           N
ATOM   3141  CA  GLN A 609       0.386  23.365 101.921  1.00 11.78           C
ATOM   3142  C   GLN A 609      -0.982  23.484 102.567  1.00 22.08           C
ATOM   3143  O   GLN A 609      -1.106  23.410 103.788  1.00 31.83           O
ATOM   3144  CB  GLN A 609       0.958  21.977 102.169  1.00 16.16           C
ATOM   3145  CG  GLN A 609       2.336  21.803 101.582  1.00 22.82           C
ATOM   3146  CD  GLN A 609       2.626  20.368 101.213  1.00 26.44           C
ATOM   3147  OE1 GLN A 609       2.506  19.466 102.043  1.00 30.04           O
ATOM   3148  NE2 GLN A 609       3.012  20.146  99.960  1.00 34.37           N
ATOM   4177  N   LEU A 737      10.811  16.643 100.700  1.00 31.38           N
ATOM   4178  CA  LEU A 737      10.324  18.015 100.557  1.00 26.08           C
ATOM   4179  C   LEU A 737      10.498  18.580  99.147  1.00 32.98           C
ATOM   4180  O   LEU A 737      11.322  18.031  98.384  1.00 22.69           O
ATOM   4181  CB  LEU A 737      11.026  18.922 101.578  1.00 28.03           C
ATOM   4182  CG  LEU A 737      10.704  20.422 101.613  1.00 21.02           C
ATOM   4183  CD1 LEU A 737      11.032  20.972 102.996  1.00 19.08           C
ATOM   4184  CD2 LEU A 737      11.488  21.164 100.534  1.00 18.74           C
ATOM   4185  OXT LEU A 737       9.810  19.576  98.824  1.00 24.22           O
END
"""

def locked(pdb, flip_state):
  """getFlippedHistidines() and the HisPlace pose for a His locked by flip_state"""
  from mmtbx.hydrogens import reduce_hydrogen, default_probe_phil
  from mmtbx.reduce import Optimizers
  m = mmtbx.model.manager(model_input=iotbx.pdb.input(source_info=None, lines=pdb),
    log=null_out())
  p = reduce_hydrogen.place_hydrogens(model=m)
  p.run()
  opt = Optimizers.Optimizer(default_probe_phil(), True, p.get_model(),
    flipStates=flip_state, fillAtomDump=False)
  poses = [l for l in opt.getInfo().splitlines() if "HisPlace at" in l and " pose " in l]
  assert len(poses) == 1, poses
  return opt.getFlippedHistidines(), poses[0]

def test_locked_unflipped_hd1_only():
  flipped, pose = locked(pdb_6vw1_his505, "1 . B HIS 505 Unflipped")
  assert "pose Unflipped HD1Placed HE2NotPlaced" in pose, pose
  assert len(flipped) == 0, pose

def test_locked_flipped_hd1_only():
  flipped, pose = locked(pdb_3j4p_his428, "1 . A HIS 428 Flipped AnglesAdjusted")
  assert "pose Flipped HD1Placed HE2NotPlaced" in pose, pose
  assert len(flipped) == 1, pose

if __name__ == "__main__":
  test_his181_keeps_ring_h()
  test_locked_unflipped_hd1_only()
  test_locked_flipped_hd1_only()
  print("OK")
