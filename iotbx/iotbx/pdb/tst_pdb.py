import iotbx.pdb.xray_structure
from cctbx import sgtbx
from cctbx.development import random_structure

def exercise_xray_structure_as_pdb_file():
  structure = random_structure.xray_structure(
    space_group_info=sgtbx.space_group_info("P 31"),
    elements=["N","C","C","O"]*2,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=00000,
    random_u_iso=0001)
  print structure.as_pdb_file(remark="Title", remarks=["Any", "Thing"]),
  print structure.as_pdb_file(remark="Title", remarks=["Any", "Thing"],
                              fractional_coordinates=0001),

def run():
  exercise_xray_structure_as_pdb_file()

if (__name__ == "__main__"):
  run()
