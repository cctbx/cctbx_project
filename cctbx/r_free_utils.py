from __future__ import absolute_import, division, print_function
from libtbx.math_utils import iceil
from libtbx.utils import null_out, Sorry
from itertools import count
import random
import sys
from six.moves import range

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

def assign_random_r_free_flags(n_refl, fraction_free, format="cns"):
  assert (fraction_free > 0) and (fraction_free < 0.5)
  from scitbx.array_family import flex
  from libtbx.math_utils import iround
  group_size = 1/(fraction_free)
  assert group_size >= 2
  if (format == "cns") or (format == "shelx"):
    result = flex.bool(n_refl, False)
    i_start = 0
    for i_group in count(1):
      i_end = min(n_refl, iround(i_group*group_size) )
      if (i_start == i_end):
        break
      if (i_end + 1 == n_refl):
        i_end += 1
      # assert i_end - i_start >= 2
      if i_end - i_start >= 2:
        result[random.randrange(i_start, i_end)] = True
      i_start = i_end
    if (format == "shelx"):
      result_ = flex.int(n_refl, 1)
      result_.set_selected(result, -1)
      result = result_
  elif (format == "ccp4"):
    result = flex.int()
    flag_max = iround(group_size) - 1
    for i in range(n_refl):
      result.append(random.randint(0, flag_max))
  return result

def assign_r_free_flags_by_shells(n_refl, fraction_free, n_bins):
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
  for i_bin in range(n_bins):
    flags.resize(min(n_refl, flags.size()+half_n_work_per_bin), False)
    flags.resize(min(n_refl, flags.size()+n_free_per_bin), True)
    flags.resize(min(n_refl, flags.size()+half_n_work_per_bin), False)
  flags.resize(n_refl, False)
  return flags

def export_r_free_flags_for_ccp4(flags, test_flag_value):
  assert (test_flag_value == True) or isinstance(test_flag_value, int)
  from scitbx.array_family import flex
  if (isinstance(flags, flex.bool)):
    test_flag_value = True
  else :
    assert isinstance(flags, flex.int)
  unique_values = set(flags)
  if len(unique_values) > 2 : # XXX: is this safe?
    return flags
  new_flags = flex.int(flags.size())
  n_free = flags.count(test_flag_value)
  if (n_free > 0):
    n_bins = iceil(flags.size() / n_free)
  else :
    n_bins = 1 # XXX dangerous!  but necessary for tiny sets
  for i in range(flags.size()):
    if flags[i] == test_flag_value :
      new_flags[i] = 0
    else :
      new_flags[i] = iceil(random.random() * (n_bins - 1))
  return new_flags

def export_r_free_flags_for_shelx(flags, test_flag_value):
  assert (test_flag_value == True) or isinstance(test_flag_value, int)
  from scitbx.array_family import flex
  if (isinstance(flags, flex.bool)):
    test_flag_value = True
  else :
    assert isinstance(flags, flex.int)
  new_flags = flex.int(flags.size(), 1)
  for i in range(flags.size()):
    if (flags[i] == test_flag_value):
      new_flags[i] = -1
  return new_flags

def looks_like_ccp4_flags(flags):
  from scitbx.array_family import flex
  assert isinstance(flags.data(), flex.int)
  return (flex.max(flags.data()) >= 4)

def adjust_fraction(miller_array, fraction, log=None):
  """
  Expand or shrink an existing set of R-free flags to match the target
  fraction.
  """
  from scitbx.array_family import flex
  assert (isinstance(miller_array.data(), flex.bool))
  assert (fraction > 0) and (fraction < 1)
  if (log is None) : log = sys.stdout
  n_refl = miller_array.data().size()
  n_free = miller_array.data().count(True)
  assert (n_refl > 0)
  current_fraction = n_free / n_refl
  print("  current fraction free: %g" % current_fraction, file=log)
  print("                 target: %g" % fraction, file=log)
  new_flags = None
  # XXX this should probably be more approximate
  if (current_fraction == fraction):
    new_flags = miller_array
  elif (current_fraction > fraction):
    # move existing flags to work set
    sub_fraction = fraction / current_fraction
    print("  shrinking old test set by a factor of %g" % sub_fraction, file=log)
    work_set = miller_array.select(~miller_array.data())
    free_set = miller_array.select(miller_array.data())
    # XXX the code in cctbx.miller requires a fraction between 0 and 0.5 -
    # since I'm not sure if this is an implementation detail or just something
    # to prevent users from doing something stupid, I'm using this clumsy
    # workaround.
    assert (sub_fraction > 0) and (sub_fraction < 1.0)
    if (sub_fraction == 0.5):
      sub_fraction = 0.501
    if (sub_fraction >= 0.5):
      free_set_new = free_set.generate_r_free_flags(
        fraction=1.0-sub_fraction,
        max_free=None,
        use_lattice_symmetry=True)
      free_set_new = free_set_new.customized_copy(
        data=~free_set_new.data())
    else :
      free_set_new = free_set.generate_r_free_flags(
        fraction=sub_fraction,
        max_free=None,
        use_lattice_symmetry=True)
    new_flags = work_set.complete_with(other=free_set_new)
  else :
    # generate additional flags
    fraction_new = (fraction - current_fraction) / (1 - current_fraction)
    print("  flagging %g of current work set for R-free" % fraction_new, file=log)
    assert (fraction_new > 0)
    work_set = miller_array.select(~miller_array.data())
    free_set = miller_array.select(miller_array.data())
    assert (work_set.indices().size() > 0)
    flags_new = work_set.generate_r_free_flags(
      fraction=fraction_new,
      max_free=None,
      use_lattice_symmetry=True)
    new_flags = flags_new.complete_with(other=free_set)#miller_array)
  assert (new_flags.indices().size() == miller_array.indices().size())
  print("    old flags: %d free (out of %d) " % (n_free, n_refl), file=log)
  print("    new flags: %d free" % new_flags.data().count(True), file=log)
  return new_flags

# XXX if the flags don't actually need extending, the original
# values will be preserved.  I think this is a good thing, but does
# this inconsistency cause problems elsewhere?
def extend_flags(
    r_free_flags,
    test_flag_value,
    array_label,
    complete_set=None,
    accumulation_callback=None,
    preserve_input_values=False,
    allow_uniform_flags=False,
    d_max=None,
    d_min=None,
    log=None):
  from scitbx.array_family import flex
  if (log is None) : log = null_out()
  assert (test_flag_value is not None)
  r_free_as_bool = get_r_free_as_bool(r_free_flags,
    test_flag_value).data()
  assert isinstance(r_free_as_bool, flex.bool)
  if(r_free_as_bool.size() == 0):
    msg = """\
WARNING: array of R-free flags (%s) is empty.""" % array_label
    print(msg, file=log)
    return r_free_flags
  fraction_free = r_free_as_bool.count(True) / r_free_as_bool.size()
  print("%s: fraction_free=%.3f" %(array_label, fraction_free), file=log)
  if (fraction_free == 0):
    if (allow_uniform_flags):
      msg = """\
WARNING: R-free flags in %s do not appear to contain a valid test, so they \
can't be extended to higher resolution.""" % array_label
      print(msg, file=log)
      return r_free_flags
    else :
      raise Sorry(("Can't extend R-free flags in %s to higher resolution "+
        "because no valid test set with flag value '%s' was found.") %
        (array_label, test_flag_value))
  if (complete_set is not None):
    missing_set = complete_set.lone_set(r_free_flags)
  else :
    tmp_d_max, tmp_d_min = r_free_flags.d_max_min()
    if (d_max is None) : d_max = tmp_d_max
    if (d_min is None) : d_min = tmp_d_min
    missing_set = r_free_flags.complete_set(d_min=d_min,
      d_max=d_max).lone_set(r_free_flags.map_to_asu())
  n_missing = missing_set.indices().size()
  print("%s: missing %d reflections" % (array_label, n_missing), file=log)
  output_array = r_free_flags
  if (n_missing != 0):
    if (n_missing <= 20):
      # FIXME: MASSIVE CHEAT necessary for tiny sets
      missing_flags = missing_set.array(data=flex.bool(n_missing,False))
    else :
      if accumulation_callback is not None :
        if not accumulation_callback(miller_array=r_free_flags,
                                     test_flag_value=test_flag_value,
                                     n_missing=n_missing,
                                     column_label=array_label):
          return r_free_flags
      missing_flags = missing_set.generate_r_free_flags(
        fraction=fraction_free,
        max_free=None,
        use_lattice_symmetry=True)
    if (preserve_input_values):
      if (looks_like_ccp4_flags(r_free_flags)):
        print("Exporting missing flags to CCP4 convention", file=log)
        exported_flags = export_r_free_flags_for_ccp4(
          flags=missing_flags.data(),
          test_flag_value=True) #test_flag_value)
        output_array = r_free_flags.concatenate(
          other=missing_flags.customized_copy(data=exported_flags))
      else :
        # XXX this is gross too - what conventions (if any) should be
        # followed here?
        work_flag_value = None
        if (test_flag_value in [1,-1]):
          work_flag_value = 0
        elif (test_flag_value == 0):
          work_flag_value = 1
        if (work_flag_value is None):
          raise Sorry(("PHENIX doesn't know how to deal with the "+
            "R-free flag convention in %s; you will need to "+
            "disable either extending the flags or preserving the "+
            "input values.") % (array_label))
        exported_flags = flex.int()
        new_flags = missing_flags.data()
        for i_seq in range(missing_flags.data().size()):
          if (new_flags[i_seq]):
            exported_flags.append(test_flag_value)
          else :
            exported_flags.append(work_flag_value)
        output_array = r_free_flags.concatenate(
          other=missing_flags.customized_copy(data=exported_flags))
    else :
      output_array = r_free_flags.concatenate(other=missing_flags)
  return output_array

def get_r_free_stats(miller_array, test_flag_value):
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
    if x > (n_ref_last + 1):
      n_bins += 1
    n_ref_last = x
  return (n_bins, n_free, sse, accu)

def get_r_free_as_bool(miller_array, test_flag_value=0):
  if miller_array.is_bool_array():
    return miller_array
  else :
    from scitbx.array_family import flex
    assert isinstance(test_flag_value, int)
    assert miller_array.is_integer_array()
    new_data = miller_array.data() == test_flag_value
    assert isinstance(new_data, flex.bool)
    return miller_array.customized_copy(data=new_data, sigmas=None)

def remediate_mismatches(array, verbose=False, log=None):
  """
  Given a set of R-free flags generated for anomalous data, detect any
  mismatches between Friedel/Bijvoet mates, and move reflections to the
  free set as needed to ensure consistency.
  """
  if (log is None) : log = null_out()
  if (not array.anomalous_flag()):
    return array
  assert array.is_bool_array()
  array = array.map_to_asu()
  orig_array = array.deep_copy()
  array = array.as_non_anomalous_array()
  merge = array.merge_equivalents(incompatible_flags_replacement=True)
  if (merge.n_incompatible_flags > 0):
    print("  %d reflections moved to test set" % \
      merge.n_incompatible_flags, file=log)
  new_flags = merge.array().generate_bijvoet_mates().common_set(
    other=orig_array)
  assert len(new_flags.indices()) == len(orig_array.indices())
  return new_flags
