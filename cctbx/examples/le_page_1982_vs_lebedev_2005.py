"""
Plot of Le Page 1982 deltas vs. Lebedev 2005 perturbations based on
random sampling of distorted unit cells compatible with the 81 2-fold
symmetry operations possible for reduced cells.
"""
from __future__ import absolute_import, division, print_function

from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
import random
from six.moves import range
from six.moves import zip

random.seed(0)

def enumerate_reduced_cell_two_folds():
  result = []
  for elements in flex.nested_loop([-1]*9,[1+1]*9):
    r = sgtbx.rot_mx(elements)
    if (r.determinant() != 1): continue
    if (r.inverse() != r): continue
    if (r.is_unit_mx()): continue
    result.append(r)
  return result

def sample(two_folds, fudge_factor, deltas, perturbations):
  for two_fold in two_folds:
    group = sgtbx.space_group()
    group.expand_smx(sgtbx.rt_mx(two_fold))
    assert group.order_z() == 2
    sym = sgtbx.space_group_info(group=group).any_compatible_crystal_symmetry(
      volume=1000)
    for i_trial in range(30):
      while True:
        uc_fudge = list(sym.unit_cell().parameters())
        for i in range(6): uc_fudge[i] *= 1+(random.random()*2-1)*fudge_factor
        try: uc_fudge = uctbx.unit_cell(uc_fudge)
        except ValueError: pass
        else: break
      deltas.append(
        two_fold.le_page_1982_delta(reduced_cell=uc_fudge, deg=True))
      perturbations.append(
        two_fold.lebedev_2005_perturbation(reduced_cell=uc_fudge))

def run():
  two_folds = enumerate_reduced_cell_two_folds()
  assert len(two_folds) == 81
  deltas = flex.double()
  perturbations = flex.double()
  for fudge_factor in [0.002, 0.01, 0.02, 0.05, 0.1]:
    sample(
      two_folds=two_folds,
      fudge_factor=fudge_factor,
      deltas=deltas,
      perturbations=perturbations)
  perm = flex.sort_permutation(data=deltas)
  deltas = deltas.select(perm)
  perturbations = perturbations.select(perm)
  f = open("le_page_1982_vs_lebedev_2005_plot", "w")
  for x,y in zip(deltas, perturbations):
    print(x, y, file=f)
  print("OK")

if (__name__ == "__main__"):
  run()
