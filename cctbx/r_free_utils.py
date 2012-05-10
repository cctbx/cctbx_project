from __future__ import division
from libtbx.math_utils import iceil
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

def assign_random_r_free_flags (n_refl, fraction_free, format="cns") :
  assert (fraction_free > 0) and (fraction_free < 0.5)
  from scitbx.array_family import flex
  from libtbx.math_utils import iround
  group_size = 1/(fraction_free)
  assert group_size >= 2
  if (format == "cns") or (format == "shelx") :
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
    if (format == "shelx") :
      result_ = flex.int(n_refl, 1)
      result_.set_selected(result, -1)
      result = result_
  elif (format == "ccp4") :
    result = flex.int()
    flag_max = iround(group_size) - 1
    for i in range(n_refl) :
      result.append(random.randint(0, flag_max))
  return result

def assign_r_free_flags_by_shells (n_refl, fraction_free, n_bins) :
  assert (fraction_free > 0) and (fraction_free < 0.5)
  assert (n_bins > 1) and (n_bins < n_refl) # XXX this should be smarter
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

def export_r_free_flags_for_ccp4 (flags, test_flag_value) :
  assert (test_flag_value == True) or isinstance(test_flag_value, int)
  from scitbx.array_family import flex
  if (isinstance(flags, flex.bool)) :
    test_flag_value = True
  else :
    assert isinstance(flags, flex.int)
  unique_values = set(flags)
  if len(unique_values) > 2 : # XXX: is this safe?
    return flags
  new_flags = flex.int(flags.size())
  n_free = flags.count(test_flag_value)
  if (n_free > 0) :
    n_bins = iceil(flags.size() / n_free)
  else :
    n_bins = 1 # XXX dangerous!  but necessary for tiny sets
  for i in range(flags.size()) :
    if flags[i] == test_flag_value :
      new_flags[i] = 0
    else :
      new_flags[i] = iceil(random.random() * (n_bins - 1))
  return new_flags

def export_r_free_flags_for_shelx (flags, test_flag_value) :
  assert (test_flag_value == True) or isinstance(test_flag_value, int)
  from scitbx.array_family import flex
  if (isinstance(flags, flex.bool)) :
    test_flag_value = True
  else :
    assert isinstance(flags, flex.int)
  new_flags = flex.int(flags.size(), 1)
  for i in range(flags.size()) :
    if (flags[i] == test_flag_value) :
      new_flags[i] = -1
  return new_flags

def looks_like_ccp4_flags (flags) :
  from scitbx.array_family import flex
  assert isinstance(flags.data(), flex.int)
  return (flex.max(flags.data()) >= 4)

def get_r_free_stats (miller_array, test_flag_value) :
  from scitbx.array_family import flex
  array = get_r_free_as_bool(miller_array, test_flag_value)
  n_free = array.data().count(True)
  accu =  array.sort(by_value="resolution").r_free_flags_accumulation()
  lr = flex.linear_regression(accu.reflection_counts.as_double(),
                              accu.free_fractions)
  assert lr.is_well_defined()
  slope = lr.slope()
  y_ideal = accu.reflection_counts.as_double() * slope
  sse = 0
  n_bins = 0
  n_ref_last = 0
  sse = flex.sum(flex.pow(y_ideal - accu.free_fractions, 2))
  for x in accu.reflection_counts :
    if x > (n_ref_last + 1) :
      n_bins += 1
    n_ref_last = x
  return (n_bins, n_free, sse, accu)

def get_r_free_as_bool (miller_array, test_flag_value=0) :
  if miller_array.is_bool_array() :
    return miller_array
  else :
    from scitbx.array_family import flex
    assert isinstance(test_flag_value, int)
    assert miller_array.is_integer_array()
    new_data = miller_array.data() == test_flag_value
    assert isinstance(new_data, flex.bool)
    return miller_array.customized_copy(data=new_data, sigmas=None)

def exercise () :
  from cctbx import miller
  from cctbx import crystal
  from cctbx import sgtbx
  from cctbx import uctbx
  from scitbx.array_family import flex
  flags_1 = assign_random_r_free_flags(n_refl=100000, fraction_free=0.05)
  assert (flags_1.count(True) == 5000)
  flags_1_ccp4 = assign_random_r_free_flags(n_refl=100000, fraction_free=0.05,
    format="ccp4")
  # XXX this is the best we can do with the current method
  assert (flags_1_ccp4.count(0) > 4000) and (flags_1_ccp4.count(0) < 6000)
  flags_1_shelx = assign_random_r_free_flags(n_refl=100000, fraction_free=0.05,
    format="shelx")
  assert (flags_1_shelx.count(-1) == 5000)
  flags_2 = assign_r_free_flags_by_shells(n_refl=100000,
    fraction_free=0.05,
    n_bins=50)
  assert (flags_2.count(True) == 5000)
  ccp4_flags = export_r_free_flags_for_ccp4(flags_1, True)
  assert (ccp4_flags.count(0) == flags_1.count(True))
  assert (flex.max(ccp4_flags) == 19)
  shelx_flags = export_r_free_flags_for_shelx(flags_1, True)
  assert ((shelx_flags==-1).all_eq((ccp4_flags==0)))
  flags_3 = assign_random_r_free_flags(n_refl=100000, fraction_free=0.025)
  assert (flags_3.count(True) == 2500)
  ccp4_flags = export_r_free_flags_for_ccp4(flags_3, True)
  assert (ccp4_flags.count(0) == flags_3.count(True))
  assert (flex.max(ccp4_flags) == 39)
  # now with an actual Miller array
  symm = crystal.symmetry(
    space_group_info=sgtbx.space_group_info("P212121"),
    unit_cell=uctbx.unit_cell((6,7,8,90,90,90)))
  set1 = miller.build_set(
    crystal_symmetry=symm,
    anomalous_flag=True,
    d_min=1.0)
  flags_4 = set1.generate_r_free_flags()
  stats = get_r_free_stats(flags_4, True)
  assert (20 <= stats[0] <= 24) # XXX is this even necessary?

if (__name__ == "__main__") :
  exercise()
  print "OK"
