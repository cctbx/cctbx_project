import iotbx.pdb.xray_structure
from cctbx.development import random_structure
from cctbx import sgtbx
from scitbx.python_utils import easy_pickle
import sys

def generate_random_f_calc(space_group_info, n_elements=10, d_min=1.5):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=False)
  structure.show_summary().show_scatterers()
  print
  print "Writing tmp.pdb"
  s = structure.as_pdb_file(
    remark="random structure",
    res_name="RND")
  open("tmp.pdb", "w").write(s)
  print
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  f_calc.show_summary()
  print
  print "Writing: tmp.phs"
  f_calc.as_phases_phs(out=open("tmp.phs", "w"))
  print

def run():
  if (len(sys.argv) != 2):
    print __doc__
    return
  generate_random_f_calc(sgtbx.space_group_info(sys.argv[1]))

if (__name__ == "__main__"):
  run()
