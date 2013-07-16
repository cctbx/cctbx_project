
from __future__ import division
from libtbx.utils import null_out
from libtbx.test_utils import show_diff
from cStringIO import StringIO

def exercise_metal_coordination () :
  pdb_in = """\
HETATM    6  PA  ANP A   1      -2.347   2.100 -17.803  1.00 10.27           P
HETATM    7  PB  ANP A   1       0.438   2.937 -18.222  1.00 10.85           P
HETATM    8  PG  ANP A   1       1.125   2.500 -21.112  1.00 10.56           P
HETATM    9  C5' ANP A   1      -1.973  -0.334 -17.074  1.00 10.24           C
HETATM   10  O5' ANP A   1      -2.737   0.857 -16.887  1.00  9.76           O
HETATM   11  C4' ANP A   1      -2.889  -1.539 -17.210  1.00 10.40           C
HETATM   12  O4' ANP A   1      -3.627  -1.699 -15.990  1.00 10.47           O
HETATM   13  C3' ANP A   1      -3.936  -1.477 -18.319  1.00  9.24           C
HETATM   14  O3' ANP A   1      -3.371  -1.834 -19.590  1.00 10.01           O
HETATM   15  C2' ANP A   1      -4.953  -2.483 -17.795  1.00  9.88           C
HETATM   16  O2' ANP A   1      -4.490  -3.826 -17.960  1.00 11.25           O
HETATM   17  C1' ANP A   1      -4.915  -2.225 -16.303  1.00 11.21           C
HETATM   18  N1  ANP A   1      -9.509  -1.559 -14.126  1.00 11.19           N
HETATM   19  O1A ANP A   1      -3.108   3.253 -17.261  1.00 10.62           O
HETATM   20  O1B ANP A   1       0.034   4.332 -18.573  1.00 10.23           O
HETATM   21  O1G ANP A   1       2.434   1.809 -21.423  1.00 12.01           O
HETATM   22  C2  ANP A   1      -8.855  -2.697 -14.422  1.00 11.32           C
HETATM   23  O2A ANP A   1      -2.510   1.726 -19.235  1.00  9.76           O
HETATM   24  O2B ANP A   1       1.600   2.764 -17.285  1.00 12.69           O
HETATM   25  O2G ANP A   1      -0.014   2.061 -22.037  1.00 10.31           O
HETATM   26  N3  ANP A   1      -7.636  -2.728 -14.988  1.00 10.83           N
HETATM   27  N3B ANP A   1       0.614   2.050 -19.591  1.00 13.13           N
HETATM   28  O3A ANP A   1      -0.783   2.257 -17.423  1.00 11.18           O
HETATM   29  O3G ANP A   1       1.269   4.014 -21.065  1.00 10.83           O
HETATM   30  C4  ANP A   1      -7.053  -1.540 -15.290  1.00 10.41           C
HETATM   31  C5  ANP A   1      -7.703  -0.266 -15.019  1.00 10.36           C
HETATM   32  C6  ANP A   1      -9.013  -0.333 -14.380  1.00 10.39           C
HETATM   33  N6  ANP A   1      -9.709   0.784 -14.079  1.00 10.34           N
HETATM   34  N7  ANP A   1      -6.889   0.732 -15.437  1.00 10.05           N
HETATM   35  C8  ANP A   1      -5.791   0.126 -15.942  1.00 10.13           C
HETATM   36  N9  ANP A   1      -5.898  -1.223 -15.869  1.00 10.39           N
HETATM   37 MG   MG  A   2      -1.727   1.886 -21.040  1.00  8.15          MG
HETATM   38 MG   MG  A   3      -0.077   5.289 -20.356  1.00  9.42          MG
HETATM    1  O   HOH A   4      -3.435   1.768 -22.164  1.00  8.59           O
HETATM    3  O   HOH A   5      -2.158   5.884 -19.807  1.00  9.00           O
HETATM    4  O   HOH A   6      -1.791   4.005 -20.881  1.00  7.86           O
HETATM   39  O   HOH A   7      -1.518  -0.155 -21.151  1.00 13.31           O
HETATM   40  O   HOH A   8      -0.256   6.271 -22.179  1.00 12.97           O
HETATM   41  O   HOH A   9       0.913   6.863 -19.633  1.00 13.62           O  """
  open("tst_geo_min_metal_coord.pdb", "w").write(pdb_in)
  params1 = """\
geometry_restraints.edits {
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    2
    atom_selection_2 = name  O2A and chain A and resname ANP and resseq    1
    distance_ideal = 2.090000
    sigma = 0.250
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    2
    atom_selection_2 = name  O2G and chain A and resname ANP and resseq    1
    distance_ideal = 2.090000
    sigma = 0.250
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    2
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    4
    distance_ideal = 2.070000
    sigma = 0.050
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    2
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    7
    distance_ideal = 2.070000
    sigma = 0.050
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    2
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    6
    distance_ideal = 2.070000
    sigma = 0.050
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    3
    atom_selection_2 = name  O3G and chain A and resname ANP and resseq    1
    distance_ideal = 2.090000
    sigma = 0.250
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    3
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    9
    distance_ideal = 2.070000
    sigma = 0.050
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    3
    atom_selection_2 = name  O1B and chain A and resname ANP and resseq    1
    distance_ideal = 2.090000
    sigma = 0.250
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    3
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    8
    distance_ideal = 2.070000
    sigma = 0.050
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    3
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    6
    distance_ideal = 2.070000
    sigma = 0.050
  }
  bond {
    action = *add
    atom_selection_1 = name MG   and chain A and resname MG and resseq    3
    atom_selection_2 = name  O   and chain A and resname HOH and resseq    5
    distance_ideal = 2.070000
    sigma = 0.050
  }
}"""
  params2 = "refinement {\n%s\n}" % params1
  open("tst_geo_min_metal_coord_1.phil", "w").write(params1)
  open("tst_geo_min_metal_coord_2.phil", "w").write(params2)
  # phenix.geometry_minimization
  from mmtbx.command_line import geometry_minimization
  args = [
    "file_name=tst_geo_min_metal_coord.pdb",
    "output_file_name_prefix=tst_geo_min_metal_coord_1",
    "write_geo_file=True",
    "tst_geo_min_metal_coord_1.phil",
  ]
  geometry_minimization.run(args=args, log=null_out())
  geo_file_1 = open("tst_geo_min_metal_coord_1.geo").read()
  args = [
    "file_name=tst_geo_min_metal_coord.pdb",
    "output_file_name_prefix=tst_geo_min_metal_coord_2",
    "write_geo_file=True",
    "tst_geo_min_metal_coord_2.phil",
  ]
  geometry_minimization.run(args=args, log=null_out())
  geo_file_2 = open("tst_geo_min_metal_coord_2.geo").read()
  show_diff(geo_file_1, geo_file_2)
  # phenix.pdbtools
  from mmtbx import pdbtools
  args = [
    "tst_geo_min_metal_coord.pdb",
    "tst_geo_min_metal_coord_1.phil",
    "model_statistics=True",
  ]
  out1 = StringIO()
  pdbtools.run(args=args, out=out1, replace_stderr=False)
  args = [
    "tst_geo_min_metal_coord.pdb",
    "tst_geo_min_metal_coord_2.phil",
    "model_statistics=True",
  ]
  out2 = StringIO()
  pdbtools.run(args=args, out=out2, replace_stderr=False)
  assert ("""      atom 1: \"HETATM   37 MG   MG  A   2 .*.    Mg  \"""" in
    out1.getvalue())
  assert ("""      atom 1: \"HETATM   37 MG   MG  A   2 .*.    Mg  \"""" in
    out2.getvalue())
  for line1, line2 in zip(out1.getvalue().splitlines(),
                          out2.getvalue().splitlines()) :
    line1 = line1.strip()
    line2 = line2.strip()
    if (line1.startswith("Date") or line1.startswith("PID") or
        line1.startswith("Command line") or line1.startswith("Time building")):
      continue
    else :
      assert not show_diff(line1, line2)

def exercise_secondary_structure () :
  pdb_in = """\
ATOM      1  N   SER A   1       0.043  -0.056   0.000  1.00  0.00
ATOM      2  CA  SER A   1       1.502  -0.020   0.000  1.00  0.00
ATOM      3  C   SER A   1       2.009   1.402   0.000  1.00  0.00
ATOM      4  O   SER A   1       1.244   2.367   0.000  1.00  0.00
ATOM      5  CB  SER A   1       2.073  -0.829  -1.192  1.00  0.00
ATOM      6  N   ASP A   2       3.321   1.549   0.000  1.00  0.00
ATOM      7  CA  ASP A   2       3.951   2.866   0.000  1.00  0.00
ATOM      8  C   ASP A   2       3.648   3.611  -1.278  1.00  0.00
ATOM      9  O   ASP A   2       3.392   4.855  -1.287  1.00  0.00
ATOM     10  CB  ASP A   2       5.477   2.730   0.242  1.00  0.00
ATOM     11  N   PRO A   3       3.664   2.893  -2.385  1.00  0.00
ATOM     12  CA  PRO A   3       3.393   3.399  -3.800  1.00  0.00
ATOM     13  C   PRO A   3       1.922   3.678  -4.096  1.00  0.00
ATOM     14  O   PRO A   3       1.642   4.698  -4.776  1.00  0.00
ATOM     15  CB  PRO A   3       3.863   2.300  -4.742  1.00  0.00
ATOM     16  N   ALA A   4       0.973   2.874  -3.653  1.00  0.00
ATOM     17  CA  ALA A   4      -0.440   3.120  -3.924  1.00  0.00
ATOM     18  C   ALA A   4      -0.904   4.394  -3.258  1.00  0.00
ATOM     19  O   ALA A   4      -1.686   5.175  -3.813  1.00  0.00
ATOM     20  CB  ALA A   4      -1.226   1.883  -3.459  1.00  0.00
ATOM     21  N   ALA A   5      -0.436   4.616  -2.044  1.00  0.00
ATOM     22  CA  ALA A   5      -0.800   5.806  -1.281  1.00  0.00
ATOM     23  C   ALA A   5      -0.280   7.057  -1.948  1.00  0.00
ATOM     24  O   ALA A   5      -0.943   8.100  -1.993  1.00  0.00
ATOM     25  CB  ALA A   5      -0.262   5.624   0.149  1.00  0.00
ATOM     26  N   LEU A   6       0.930   6.974  -2.469  1.00  0.00
ATOM     27  CA  LEU A   6       1.561   8.105  -3.144  1.00  0.00
ATOM     28  C   LEU A   6       0.807   8.474  -4.399  1.00  0.00
ATOM     29  O   LEU A   6       0.621   9.654  -4.732  1.00  0.00
ATOM     30  CB  LEU A   6       3.050   7.789  -3.463  1.00  0.00
ATOM     31  N   LYS A   7       0.357   7.464  -5.122  1.00  0.00
ATOM     32  CA  LYS A   7      -0.388   7.671  -6.360  1.00  0.00
ATOM     33  C   LYS A   7      -1.709   8.352  -6.090  1.00  0.00
ATOM     34  O   LYS A   7      -2.158   9.234  -6.836  1.00  0.00
ATOM     35  CB  LYS A   7      -0.602   6.306  -7.069  1.00  0.00
ATOM     36  N   ARG A   8      -2.356   7.951  -5.011  1.00  0.00
ATOM     37  CA  ARG A   8      -3.644   8.521  -4.625  1.00  0.00
ATOM     38  C   ARG A   8      -3.504   9.981  -4.267  1.00  0.00
ATOM     39  O   ARG A   8      -4.355  10.812  -4.578  1.00  0.00
ATOM     40  CB  ARG A   8      -4.223   7.699  -3.441  1.00  0.00
ATOM     41  N   ALA A   9      -2.427  10.307  -3.578  1.00  0.00
ATOM     42  CA  ALA A   9      -2.161  11.682  -3.163  1.00  0.00
ATOM     43  C   ALA A   9      -1.934  12.574  -4.359  1.00  0.00
ATOM     44  O   ALA A   9      -2.376  13.729  -4.409  1.00  0.00
ATOM     45  CB  ALA A   9      -0.959  11.654  -2.202  1.00  0.00
ATOM     46  N   ARG A  10      -1.225  12.054  -5.343  1.00  0.00
ATOM     47  CA  ARG A  10      -0.927  12.801  -6.561  1.00  0.00
ATOM     48  C   ARG A  10      -2.188  13.100  -7.335  1.00  0.00
ATOM     49  O   ARG A  10      -2.361  14.174  -7.908  1.00  0.00
ATOM     50  CB  ARG A  10       0.084  11.988  -7.416  1.00  0.00
ATOM     51  N   ASN A  11      -3.083  12.131  -7.381  1.00  0.00
ATOM     52  CA  ASN A  11      -4.349  12.279  -8.093  1.00  0.00
ATOM     53  C   ASN A  11      -5.209  13.341  -7.451  1.00  0.00
ATOM     54  O   ASN A  11      -5.871  14.135  -8.118  1.00  0.00
ATOM     55  CB  ASN A  11      -5.091  10.911  -8.123  1.00  0.00
"""
  open("tst_geo_min_ss_phil.pdb", "w").write(pdb_in)
  from mmtbx.command_line import secondary_structure_restraints
  out = StringIO()
  args = [
    "tst_geo_min_ss_phil.pdb",
    "use_ksdssp=False",
  ]
  secondary_structure_restraints.run(args, out=out, log=null_out())
  # phenix.geometry_minimization
  from mmtbx.command_line import geometry_minimization
  open("tst_geo_min_ss_phil.phil", "w").write(out.getvalue())
  args = [
    "tst_geo_min_ss_phil.pdb",
    "tst_geo_min_ss_phil.phil",
    "secondary_structure_restraints=True",
    "secondary_structure.input.find_automatically=False",
  ]
  geometry_minimization.run(args=args, log=null_out())
  # TODO phenix.pdbtools?

if (__name__ == "__main__") :
  exercise_metal_coordination()
  exercise_secondary_structure()
  print "OK"
