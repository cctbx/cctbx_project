from cctbx import sgtbx

class group_search(object):

  def __init__(self, modulus=2):
    self.ext = sgtbx.lattice_symmetry_group_search(modulus)

  def n_potential_axes(self):
    return self.ext.n_potential_axes()

  def __call__(self, minimum_cell, max_delta=3.,
               only_test_generators=True):
    return self.ext(minimum_cell, max_delta, only_test_generators)

group = group_search()

def find_max_delta(minimum_cell, group, modulus=2):
  return sgtbx.lattice_symmetry_find_max_delta(minimum_cell, group, modulus)
