from __future__ import absolute_import, division, print_function
import random
from cctbx import crystal
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
from mmtbx.refinement.real_space import fit_residues
import mmtbx.utils

class symmetry_holder(object):
  # what fit_residues.run._neighbors_of uses of run
  def __init__(self, cs):
    self.crystal_symmetry = cs
    self.special_position_settings = crystal.special_position_settings(
      crystal_symmetry=cs)

def all_atoms(h, sites_cart, i_selection, radius):
  return list(h.special_position_settings.pair_generator(
    sites_cart=sites_cart, distance_cutoff=radius).neighbors_of(
      primary_selection=flex.bool(sites_cart.size(), i_selection)).iselection())

def neighbors_of(h, sites_cart, i_selection, radius):
  return list(fit_residues.run._neighbors_of(h, sites_cart=sites_cart,
    i_selection=i_selection, radius=radius))

def exercise_neighbors_of():
  # the restricted search gives pair_generator's result for all atoms,
  # also next to atoms that pair_generator moves onto a symmetry element
  # (here a 2-fold along y through x = z = 0), and for an empty selection
  h = symmetry_holder(crystal.symmetry(unit_cell=(30, 30, 30, 90, 90, 90),
    space_group_symbol="P2"))
  for sites, isel, radius in [
      ([(0.2, 10, 0), (0, 10, 1.49)], [1], 1.5),
      ([(0.2, 10, 0), (0, 10, 1.49), (0, 10, 3.0), (0, 10, 4.51)], [1, 2, 3], 1.5),
      ([(0.0, 10, 1.6), (0.2, 10, 0.0)], [1], 1.5),
      ([(0.2, 10, 0), (0, 10, 4.99)], [1], 5.0),
      ([(0.3, 10, 0), (0, 10, 1.48), (0, 10, 6.0)], [1, 2], 1.5),
      ([(1, 2, 3), (4, 5, 6)], [], 3.0)]:
    sites, isel = flex.vec3_double(sites), flex.size_t(isel)
    assert neighbors_of(h, sites, isel, radius) == \
      all_atoms(h, sites, isel, radius), (list(sites), list(isel))
  # general positions, atoms anywhere (also outside the cell)
  rnd = random.Random(0)
  for symbol, cell in [("P 1", (31, 42, 53, 75, 95, 110)),
      ("P 1 21 1", (35, 44, 51, 90, 105, 90)), ("C 1 2 1", (60, 40, 35, 90, 115, 90)),
      ("P 61 2 2", (40, 40, 90, 90, 90, 120)), ("R 3 :H", (50, 50, 70, 90, 90, 120)),
      ("F 4 3 2", (45, 45, 45, 90, 90, 90))]:
    h = symmetry_holder(crystal.symmetry(unit_cell=cell,
      space_group_symbol=symbol))
    for trial in range(10):
      sites = flex.vec3_double([(rnd.uniform(-60, 90), rnd.uniform(-60, 90),
        rnd.uniform(-60, 90)) for i in range(150)])
      centre = sites[0]
      isel = flex.size_t([0] + [i for i in range(1, 150)
        if sum((a - b)**2 for a, b in zip(sites[i], centre)) < 9][:5])
      for radius in (1.5, 5.0):
        assert neighbors_of(h, sites, isel, radius) == \
          all_atoms(h, sites, isel, radius), (symbol, trial)

def exercise_approx_equal_fast():
  # the same pass/fail as approx_equal, without output
  from six.moves import cStringIO as StringIO
  rnd = random.Random(1)
  specials = [0., -0., 1.e300, -1.e300, float("inf"), float("-inf"),
    float("nan"), 5.e-324, 1.e-6]
  for trial in range(300):
    n = rnd.randint(0, 6)
    a = [rnd.choice(specials) if rnd.random() < 0.3 else rnd.uniform(-2, 2)
         for i in range(n)]
    b = [x + rnd.choice([0, 1.e-7, 1.e-6, 1.e-5, 1.]) if rnd.random() < 0.7
         else rnd.choice(specials) for x in a]
    for eps in (0, 1.e-6, 1.e-3):
      expected = approx_equal(flex.double(a), flex.double(b), eps=eps,
        out=StringIO())
      assert mmtbx.utils._approx_equal_fast(flex.double(a), flex.double(b),
        eps) == bool(expected), (a, b, eps)

def run():
  exercise_neighbors_of()
  exercise_approx_equal_fast()
  print("OK")

if (__name__ == "__main__"):
  run()
