
from __future__ import absolute_import, division, print_function
from libtbx.utils import null_out
from libtbx.test_utils import show_diff
from six.moves import cStringIO as StringIO
from six.moves import zip

def exercise_metal_coordination():
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
                          out2.getvalue().splitlines()):
    line1 = line1.strip()
    line2 = line2.strip()
    if (line1.startswith("Date") or line1.startswith("PID") or
        line1.startswith("Command line") or line1.startswith("Time")):
      continue
    else :
      assert not show_diff(line1, line2)

def exercise_secondary_structure():
  # what is tested here? Not a single assert or check.
  pdb_in_p1 = """\
CRYST1   21.348   25.003   18.365  90.00  90.00  90.00 P 1
ATOM      1  N   SER A   1      10.914   5.773  13.123  1.00  0.00
ATOM      2  CA  SER A   1      12.373   5.809  13.123  1.00  0.00
ATOM      3  C   SER A   1      12.880   7.231  13.123  1.00  0.00
ATOM      4  O   SER A   1      12.115   8.196  13.123  1.00  0.00
ATOM      5  CB  SER A   1      12.944   5.000  11.931  1.00  0.00
ATOM      6  N   ASP A   2      14.192   7.378  13.123  1.00  0.00
ATOM      7  CA  ASP A   2      14.822   8.695  13.123  1.00  0.00
ATOM      8  C   ASP A   2      14.519   9.440  11.845  1.00  0.00
ATOM      9  O   ASP A   2      14.263  10.684  11.836  1.00  0.00
ATOM     10  CB  ASP A   2      16.348   8.559  13.365  1.00  0.00
ATOM     11  N   PRO A   3      14.535   8.722  10.738  1.00  0.00
ATOM     12  CA  PRO A   3      14.264   9.228   9.323  1.00  0.00
ATOM     13  C   PRO A   3      12.793   9.507   9.027  1.00  0.00
ATOM     14  O   PRO A   3      12.513  10.527   8.347  1.00  0.00
ATOM     15  CB  PRO A   3      14.734   8.129   8.381  1.00  0.00
ATOM     16  N   ALA A   4      11.844   8.703   9.470  1.00  0.00
ATOM     17  CA  ALA A   4      10.431   8.949   9.199  1.00  0.00
ATOM     18  C   ALA A   4       9.967  10.223   9.865  1.00  0.00
ATOM     19  O   ALA A   4       9.185  11.004   9.310  1.00  0.00
ATOM     20  CB  ALA A   4       9.645   7.712   9.664  1.00  0.00
ATOM     21  N   ALA A   5      10.435  10.445  11.079  1.00  0.00
ATOM     22  CA  ALA A   5      10.071  11.635  11.842  1.00  0.00
ATOM     23  C   ALA A   5      10.591  12.886  11.175  1.00  0.00
ATOM     24  O   ALA A   5       9.928  13.929  11.130  1.00  0.00
ATOM     25  CB  ALA A   5      10.609  11.453  13.272  1.00  0.00
ATOM     26  N   LEU A   6      11.801  12.803  10.654  1.00  0.00
ATOM     27  CA  LEU A   6      12.432  13.934   9.979  1.00  0.00
ATOM     28  C   LEU A   6      11.678  14.303   8.724  1.00  0.00
ATOM     29  O   LEU A   6      11.492  15.483   8.391  1.00  0.00
ATOM     30  CB  LEU A   6      13.921  13.618   9.660  1.00  0.00
ATOM     31  N   LYS A   7      11.228  13.293   8.001  1.00  0.00
ATOM     32  CA  LYS A   7      10.483  13.500   6.763  1.00  0.00
ATOM     33  C   LYS A   7       9.162  14.181   7.033  1.00  0.00
ATOM     34  O   LYS A   7       8.713  15.063   6.287  1.00  0.00
ATOM     35  CB  LYS A   7      10.269  12.135   6.054  1.00  0.00
ATOM     36  N   ARG A   8       8.515  13.780   8.112  1.00  0.00
ATOM     37  CA  ARG A   8       7.227  14.350   8.498  1.00  0.00
ATOM     38  C   ARG A   8       7.367  15.810   8.856  1.00  0.00
ATOM     39  O   ARG A   8       6.516  16.641   8.545  1.00  0.00
ATOM     40  CB  ARG A   8       6.648  13.528   9.682  1.00  0.00
ATOM     41  N   ALA A   9       8.444  16.136   9.545  1.00  0.00
ATOM     42  CA  ALA A   9       8.710  17.511   9.960  1.00  0.00
ATOM     43  C   ALA A   9       8.937  18.403   8.764  1.00  0.00
ATOM     44  O   ALA A   9       8.495  19.558   8.714  1.00  0.00
ATOM     45  CB  ALA A   9       9.912  17.483  10.921  1.00  0.00
ATOM     46  N   ARG A  10       9.646  17.883   7.780  1.00  0.00
ATOM     47  CA  ARG A  10       9.944  18.630   6.562  1.00  0.00
ATOM     48  C   ARG A  10       8.683  18.929   5.788  1.00  0.00
ATOM     49  O   ARG A  10       8.510  20.003   5.215  1.00  0.00
ATOM     50  CB  ARG A  10      10.955  17.817   5.707  1.00  0.00
ATOM     51  N   ASN A  11       7.788  17.960   5.742  1.00  0.00
ATOM     52  CA  ASN A  11       6.522  18.108   5.030  1.00  0.00
ATOM     53  C   ASN A  11       5.662  19.170   5.672  1.00  0.00
ATOM     54  O   ASN A  11       5.000  19.964   5.005  1.00  0.00
ATOM     55  CB  ASN A  11       5.780  16.740   5.000  1.00  0.00
TER
END
  """
  open("tst_geo_min_ss_phil.pdb", "w").write(pdb_in_p1)
  from mmtbx.command_line import secondary_structure_restraints
  out = StringIO()
  args = [
    "tst_geo_min_ss_phil.pdb",
    # "use_ksdssp=False",
  ]
  secondary_structure_restraints.run(args, out=out, log=null_out())
  # phenix.geometry_minimization
  from mmtbx.command_line import geometry_minimization
  open("tst_geo_min_ss_phil.phil", "w").write(out.getvalue())
  args = [
    "tst_geo_min_ss_phil.pdb",
    "tst_geo_min_ss_phil.phil",
    "secondary_structure.enabled=True",
  ]
  geometry_minimization.run(args=args, log=null_out())
  # TODO phenix.pdbtools?

if (__name__ == "__main__"):
  # exercise_metal_coordination()
  exercise_secondary_structure()
  print("OK")
