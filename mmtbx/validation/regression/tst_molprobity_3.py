from __future__ import absolute_import, division, print_function
from libtbx import easy_run
import libtbx.load_env

pdb_str_1 = """
CRYST1   69.982  108.450  152.761 106.22 101.81  98.32 P 1
ATOM  21640  CA  ILE D 873      30.355 -45.714 -17.058  1.00 37.13           C
ATOM  21641  CB  ILE D 873      29.087 -46.471 -16.558  1.00 37.66           C
ATOM  21642  CG1 ILE D 873      29.392 -47.193 -15.226  1.00 39.39           C
ATOM  21643  CD1 ILE D 873      28.156 -47.608 -14.464  1.00 38.31           C
ATOM  21644  CG2 ILE D 873      28.645 -47.549 -17.518  1.00 37.63           C
ATOM  21645  C   ILE D 873      30.411 -45.513 -18.582  1.00 40.34           C
ATOM  21646  O   ILE D 873      29.507 -44.925 -19.217  1.00 37.47           O
ATOM  21647  N   VAL D 874      31.486 -46.033 -19.167  1.00 46.20           N
ATOM  21648  CA  VAL D 874      31.874 -45.719 -20.545  1.00 53.20           C
ATOM  21649  CB  VAL D 874      33.422 -45.841 -20.767  1.00 47.96           C
ATOM  21650  C   VAL D 874      31.111 -46.672 -21.465  1.00 51.88           C
ATOM  21651  O   VAL D 874      29.913 -46.492 -21.675  1.00 47.01           O
TER   21652      VAL D 874
HETATM21653  O   HOH E   1      25.038 -21.838 -79.915  1.00 30.89           O
HETATM21654  O   HOH E   2      59.047 -10.310 -77.729  1.00 32.39           O
HETATM21655  O   HOH E   3      69.611   0.997 -43.917  1.00 19.06           O
HETATM21656  O   HOH E   4      54.107 -16.143 -62.221  1.00 26.11           O
HETATM21657  O   HOH E   5      13.823 -63.943 -12.236  1.00 22.79           O
HETATM21658  O   HOH E   6      59.530 -22.000 -79.110  1.00 30.64           O
HETATM21659  O   HOH E   7      31.927 -44.110 -81.232  1.00 24.49           O
END
"""

pdb_str_2 = """
CRYST1   69.982  108.450  152.761 106.22 101.81  98.32 P 1
HETATM21659  O   HOH E   7      31.927 -44.110 -81.232  1.00 24.49           O
END
"""

pdb_str_3 = """\
CRYST1    1.000    1.000    1.000  90.00  90.00  90.00 P 1           1
ATOM      1  N   LEU B 111     -22.437  63.276  48.231  1.00 38.59           N
ATOM      2  CA  LEU B 111     -22.866  63.139  46.796  1.00 38.96           C
ATOM      3  C   LEU B 111     -22.155  62.012  46.022  1.00 39.03           C
ATOM      4  O   LEU B 111     -21.022  61.639  46.348  1.00 38.95           O
ATOM      5  CB  LEU B 111     -22.696  64.458  46.041  1.00 38.99           C
ATOM      6  CG  LEU B 111     -23.630  65.621  46.381  1.00 39.42           C
ATOM      7  CD1 LEU B 111     -23.034  66.475  47.484  1.00 39.12           C
ATOM      8  CD2 LEU B 111     -23.909  66.475  45.146  1.00 39.58           C
ATOM      9  N   TYR B 112     -22.844  61.491  44.997  1.00 39.14           N
ATOM     10  CA  TYR B 112     -22.380  60.363  44.173  1.00 39.10           C
ATOM     11  C   TYR B 112     -22.735  60.561  42.692  1.00 39.27           C
ATOM     12  O   TYR B 112     -23.901  60.786  42.361  1.00 39.37           O
ATOM     13  CB  TYR B 112     -23.002  59.035  44.649  1.00 38.91           C
ATOM     14  CG  TYR B 112     -22.527  58.549  46.004  1.00 38.64           C
ATOM     15  CD1 TYR B 112     -23.191  58.926  47.170  1.00 39.16           C
ATOM     16  CD2 TYR B 112     -21.417  57.712  46.121  1.00 38.27           C
ATOM     17  CE1 TYR B 112     -22.756  58.502  48.424  1.00 39.01           C
ATOM     18  CE2 TYR B 112     -20.978  57.269  47.370  1.00 38.32           C
ATOM     19  CZ  TYR B 112     -21.655  57.667  48.518  1.00 38.85           C
ATOM     20  OH  TYR B 112     -21.238  57.244  49.764  1.00 38.58           O
END
"""

def exercise_00(prefix="tst_molprobity_3_exercise_00"):
  for i, pdb_str in enumerate([pdb_str_1, pdb_str_2]):
    of = open("%s.pdb"%prefix, "w")
    print(pdb_str, file=of)
    of.close()
    cmd = " ".join([
      "phenix.fmodel",
      "%s.pdb"%prefix,
      "high_res=10",
      "type=real r_free=0.1 label=F-obs",
      "output.file_name=%s_%d.mtz" % (prefix, i),
      "> %s.zlog"%prefix])
    assert not easy_run.call(cmd)
    cmd = "phenix.molprobity %s.pdb %s_%d.mtz > %s.zlog"%(prefix,prefix,i,prefix)
    assert not easy_run.call(cmd)

def exercise_01(prefix="tst_molprobity_3_exercise_01"):
  of = open("%s.pdb"%prefix, "w")
  print(pdb_str_3, file=of)
  of.close()
  cmd = "phenix.molprobity %s.pdb > %s.zlog"%(prefix,prefix)
  r = easy_run.fully_buffered(cmd)
  assert r.stderr_lines[0]=="Sorry: Crystal symmetry is missing or cannot be extracted."

if (__name__ == "__main__"):
  if libtbx.env.has_module("phenix"):
    exercise_00()
    # exercise_01() disabling because such CRYST1 is working now...
    print("OK")
  else:
    print("Skipped: Requires phenix module")
