from cctbx.development import random_structure
from cctbx import miller
from cctbx.array_family import flex
import random

def random_f_calc(space_group_info, n_scatterers, d_min, anomalous_flag,
                  verbose=0):
  if (anomalous_flag and space_group_info.group().is_centric()):
    return None
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers,
    volume_per_atom=500,
    min_distance=2.,
    general_positions_only=True)
  if (0 or verbose):
    structure.show_summary().show_scatterers()
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag).f_calc()
  f_calc = miller.array(
    miller_set=f_calc,
    data=f_calc.data()/flex.mean(flex.abs(f_calc.data())))
  if (f_calc.anomalous_flag()):
    selection = flex.bool(f_calc.indices().size(), True)
    for i in xrange(f_calc.indices().size()//10):
      j = random.randrange(f_calc.indices().size())
      selection[j] = False
    f_calc = f_calc.select(selection)
  return f_calc
