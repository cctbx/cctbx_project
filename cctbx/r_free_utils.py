
from itertools import count
import random

generate_r_free_params_str = """\
    fraction = 0.1
      .type=float
      .short_caption = Fraction of reflections in test set
      .expert_level=0
    max_free = 2000
      .type=int
      .short_caption = Maximum number of reflections in test set
      .expert_level=2
    lattice_symmetry_max_delta = 5
      .type=float
      .expert_level=2
    use_lattice_symmetry = True
      .type=bool
      .short_caption = Use lattice symmetry to generate test set
      .expert_level=0
    use_dataman_shells = False
      .type = bool
      .short_caption = Assign test set in thin resolution shells
      .help = Used to avoid biasing of the test set by certain types of \
        non-crystallographic symmetry.
    n_shells = 20
      .type = int
      .short_caption = Number of resolution shells
"""

def assign_random_r_free_flags (n_refl, fraction_free) :
  from scitbx.array_family import flex
  from libtbx.math_utils import iround
  group_size = 1/(fraction_free)
  assert group_size >= 2
  result = flex.bool(n_refl, False)
  i_start = 0
  for i_group in count(1):
    i_end = min(n_refl, iround(i_group*group_size) )
    if (i_start == i_end):
      break
    if (i_end + 1 == n_refl):
      i_end += 1
    assert i_end - i_start >= 2
    result[random.randrange(i_start, i_end)] = True
    i_start = i_end
  return result

def assign_r_free_flags_by_shells (n_refl, fraction_free, n_bins) :
  from scitbx.array_family import flex
  from libtbx.math_utils import iround
  n_free = iround(n_refl * fraction_free)
  n_per_bin = iround(n_refl / n_bins)
  half_n_work_per_bin = iround((n_refl-n_free) / n_bins / 2)
  n_free_per_bin = n_per_bin - 2*half_n_work_per_bin
  flags = flex.bool()
  flags.reserve(n_refl)
  for i_bin in xrange(n_bins):
    flags.resize(min(n_refl, flags.size()+half_n_work_per_bin), False)
    flags.resize(min(n_refl, flags.size()+n_free_per_bin), True)
    flags.resize(min(n_refl, flags.size()+half_n_work_per_bin), False)
  flags.resize(n_refl, False)
  return flags
