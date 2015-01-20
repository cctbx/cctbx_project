from __future__ import division
from libtbx import easy_run

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

def exercise_00(prefix="tst_molprobity_3_exercise_00"):
  for pdb_str in [pdb_str_1, pdb_str_2]:
    of = open("%s.pdb"%prefix, "w")
    print >> of, pdb_str
    of.close()
    cmd = " ".join([
      "phenix.fmodel",
      "%s_in.pdb"%prefix,
      "high_res=10",
      "type=real r_free=0.1 label=F-obs",
      "output.file_name=%s.mtz"%prefix,
      "> %s.zlog"%prefix])
    easy_run.call(cmd)
    cmd = "phenix.molprobity %s.pdb %s.mtz > %s.zlog"%(prefix,prefix,prefix)
    easy_run.call(cmd)

if (__name__ == "__main__") :
  exercise_00()
  print "OK"
