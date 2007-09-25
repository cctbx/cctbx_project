from __future__ import division

from cctbx import sgtbx
from cctbx import miller
from cctbx import maptbx
from cctbx.development import random_structure

import smtbx.ab_initio.charge_flipping as cflip

def exercise_randomly(space_group_info, elements,
                      special_handling_of_weak_reflections,
                      anomalous_flag, d_min=0.8, grid_resolution_factor=1./2,
                      ):
  structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    use_u_iso=False,
    use_u_aniso=False,
  )
  for s in structure.scatterers():
    assert s.flags.use_u_iso() and s.u_iso == 0
    assert not s.flags.use_u_aniso()
  f_indices = miller.build_set(
    crystal_symmetry=structure,
    anomalous_flag=anomalous_flag,
    d_min=d_min)
  f_correct = f_indices.structure_factors_from_scatterers(
    xray_structure=structure,
    algorithm="direct").f_calc()
  f_obs = abs(f_correct)
  if special_handling_of_weak_reflections:
    flipped = cflip.charge_flipping_iterator_tweaked_for_weak_reflections(f_obs)
  else:
    flipped = cflip.charge_flipping_iterator(f_obs)
  print "start:"
  print "fraction of positive pixels=%.3f" % (
    (flipped.rho > 0).count(True)/flipped.rho.size())
  for i, state in enumerate(flipped):
    if i == 10: state.adjust_delta()
    r1 = state.f_obs.r1_factor(state.g, assume_index_matching=True)
    pos = (state.rho > 0).count(True) / state.rho.size()
    print "%i: delta=%.5f" % (i,state.delta)
    indent = " "*(len("%i" % i) + 1)
    print indent, "R1=%.5f" % r1
    print indent, "fraction of positive pixels=%.3f" % pos
    print indent, "total charge=%.3f" % state.g_000


def run():
  exercise_randomly(space_group_info=sgtbx.space_group_info("P31"),
                    elements=["C"]*5 + ["O"]*2 + ["N"],
                    anomalous_flag=False,
                    special_handling_of_weak_reflections=True
                  )

if __name__ == '__main__':
  run()
