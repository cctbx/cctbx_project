from __future__ import absolute_import, division, print_function
from mmtbx.programs import fmodel
from iotbx.cli_parser import run_program
from libtbx.utils import null_out, Sorry
import os

def exercise():
  model_fn = "tmp_fmodel_fake_p1.pdb"
  with open(model_fn, "w") as f:
    f.write("""\
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
  args = [model_fn, "high_resolution=2",
    "output.file_name=tmp_fmodel_fake_p1.mtz"]
  try:
    run_program(program_class=fmodel.Program, args=args, logger=null_out())
  except Sorry as s:
    assert('Symmetry information in model file is incomplete or missing'
      in str(s))
  else:
    raise Exception_expected
  args.append("generate_fake_p1_symmetry=True")
  r = run_program(program_class=fmodel.Program, args=args, logger=null_out())
  assert os.path.isfile(r.output_file)
  from iotbx import crystal_symmetry_from_any
  symm = crystal_symmetry_from_any.extract_from(r.output_file)
  assert (str(symm.space_group_info()) == "P 1")
  args.append(r.output_file)
  args.append("labels.name=FMODEL,PHIFMODEL")
  try:
    run_program(program_class=fmodel.Program, args=args, logger=null_out())
  except Sorry as s:
    assert('high_resolution and low_resolution must be undefined if reflection'
      in str(s))
  else:
    raise Exception_expected

  # Clean up files
  os.remove(r.output_file)
  os.remove(model_fn)

  print("OK")

if (__name__ == "__main__"):
  exercise()
