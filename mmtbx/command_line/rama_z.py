from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME mmtbx.rama_z
# LIBTBX_SET_DISPATCHER_NAME phenix.rama_z

from mmtbx.programs import rama_z
from iotbx.cli_parser import run_program

result = run_program(rama_z.Program)
if result is None:
  print("Calculation of z-score failed for some reason")
else:
  for k in ["whole", "helix", "sheet", "loop"]:
    rc = k[0].upper()
    v = result.get(rc, None)
    if v is None:
      print("z-score %-5s: None, residues: %d" % (k, result['residue_counts'][rc]))
    else:
      print("z-score %-5s: %6.3f, residues: %d" % (k, v, result['residue_counts'][rc]))
