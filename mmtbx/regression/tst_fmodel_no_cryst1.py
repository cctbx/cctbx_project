
from __future__ import division
from __future__ import print_function
from libtbx import easy_run
import os

def exercise():
  open("tmp_fmodel_fake_p1.pdb", "w").write("""\
ATOM     47  N   TYR A   7       8.292   1.817   6.147  1.00 14.70           N
ATOM     48  CA  TYR A   7       9.159   2.144   7.299  1.00 15.18           C
ATOM     49  C   TYR A   7      10.603   2.331   6.885  1.00 15.91           C
ATOM     50  O   TYR A   7      11.041   1.811   5.855  1.00 15.76           O
ATOM     51  CB  TYR A   7       9.061   1.065   8.369  1.00 15.35           C
ATOM     52  CG  TYR A   7       7.665   0.929   8.902  1.00 14.45           C
ATOM     53  CD1 TYR A   7       6.771   0.021   8.327  1.00 15.68           C
ATOM     54  CD2 TYR A   7       7.210   1.756   9.920  1.00 14.80           C
ATOM     55  CE1 TYR A   7       5.480  -0.094   8.796  1.00 13.46           C
ATOM     56  CE2 TYR A   7       5.904   1.649  10.416  1.00 14.33           C
ATOM     57  CZ  TYR A   7       5.047   0.729   9.831  1.00 15.09           C
ATOM     58  OH  TYR A   7       3.766   0.589  10.291  1.00 14.39           O
ATOM     59  OXT TYR A   7      11.358   2.999   7.612  1.00 17.49           O
""")
  args = ["phenix.fmodel", "tmp_fmodel_fake_p1.pdb", "high_resolution=2",
    "output.file_name=tmp_fmodel_fake_p1.mtz"]
  result = easy_run.fully_buffered(args)
  assert (result.return_code != 0) and (len(result.stderr_lines) > 0)
  args.append("generate_fake_p1_symmetry=True")
  result = easy_run.fully_buffered(args).raise_if_errors()
  assert (result.return_code == 0)
  assert os.path.isfile("tmp_fmodel_fake_p1.mtz")
  from iotbx import crystal_symmetry_from_any
  symm = crystal_symmetry_from_any.extract_from("tmp_fmodel_fake_p1.mtz")
  assert (str(symm.space_group_info()) == "P 1")
  # FIXME this should fail but doesn't due to a bug in the program
  #args.append("reference_file=tmp_fmodel_fake_p1.mtz")
  #args.append("data_column_label=FMODEL,PHIFMODEL")
  #result = easy_run.fully_buffered(args).raise_if_errors()
  #print result.return_code
  print("OK")

if (__name__ == "__main__"):
  exercise()
