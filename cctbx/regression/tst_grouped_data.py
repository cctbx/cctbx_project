from __future__ import absolute_import, division, print_function
from cctbx import miller
from cctbx import xray
from cctbx.development import random_structure
from cctbx.development import debug_utils
from cctbx.array_family import flex
import sys

def exercise(space_group_info, anomalous_flag,
             n_scatterers=8, d_min=2, verbose=0):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["const"]*n_scatterers)
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=anomalous_flag).f_calc()
  f = abs(f_calc)
  f = miller.array(miller_set=f, data=f.data(), sigmas=flex.sqrt(f.data()))
  f = f.f_as_f_sq()
  g = f.expand_to_p1()
  merger_p1 = xray.merger( g.indices(),
                              g.data(),
                              g.sigmas(),
                              g.space_group(),
                              g.anomalous_flag(),
                              g.unit_cell() )
  p1_bic = merger_p1.bic()
  p1_r = merger_p1.r_abs()

  merger_nat = xray.merger( g.indices(),
                              g.data(),
                              g.sigmas(),
                              f.space_group(),
                              g.anomalous_flag(),
                              g.unit_cell() )
  nat_bic = merger_nat.bic()
  nat_r = merger_nat.r_abs()
  assert nat_bic >= p1_bic
  assert p1_r <= 1e-8


def run_call_back(flags, space_group_info):
  for anomalous_flag in (False, True):
    exercise(space_group_info, anomalous_flag)

def run():
  debug_utils.parse_options_loop_space_groups(sys.argv[1:], run_call_back)
  print("OK")

if (__name__ == "__main__"):
  run()
