from cctbx import miller
from cctbx import crystal
from cctbx import sgtbx
import cctbx.sgtbx.lattice_symmetry
import cctbx.sgtbx.bravais_types
import cctbx.sgtbx.subgroups
import cctbx.sgtbx.cosets
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import random
import sys

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def exercise(space_group_info, anomalous_flag, verbose):
  crystal_symmetry_ri = crystal.symmetry(
    unit_cell=space_group_info.any_compatible_unit_cell(volume=1000),
    space_group_info=space_group_info) \
      .minimum_cell() \
      .reflection_intensity_symmetry(anomalous_flag=anomalous_flag)
  #
  crystal_symmetry_p1 = crystal_symmetry_ri.cell_equivalent_p1()
  miller_set = miller.build_set(
    crystal_symmetry=crystal_symmetry_p1,
    anomalous_flag=anomalous_flag,
    d_min=1)
  miller_array_p1 = miller.array(
    miller_set=miller_set,
    data=flex.random_double(size=miller_set.indices().size()))
  assert miller_array_p1.map_to_asu().indices() \
         .all_eq(miller_array_p1.indices())
  #
  lattice_group = sgtbx.lattice_symmetry.group(
    crystal_symmetry_ri.unit_cell(),
    max_delta=0.1)
  miller_array_subs = []
  miller_array_sub_as = []
  miller_array_sub_bs = []
  coset_decompositions = []
  subgroups = sgtbx.subgroups.anomalous_reflection_intensity_primitive_cell(
    space_group=lattice_group)
  for subgroup in subgroups:
    subgroup_info = sgtbx.space_group_info(group=subgroup)
    miller_array_sub = miller_array_p1.customized_copy(
      space_group_info=subgroup_info) \
        .merge_equivalents().array()
    miller_array_sub_a = miller_array_sub.customized_copy(
      data=flex.random_double(size=miller_array_sub.indices().size())) \
        .expand_to_p1() \
        .map_to_asu()
    miller_array_sub_b = miller_array_sub.customized_copy(
      data=flex.random_double(size=miller_array_sub.indices().size())) \
        .expand_to_p1() \
        .map_to_asu()
    miller_array_subs.append(miller_array_sub)
    miller_array_sub_as.append(miller_array_sub_a)
    miller_array_sub_bs.append(miller_array_sub_b)
    #
    for s in subgroup_info.group():
      cb_op = sgtbx.change_of_basis_op(s)
      cb = miller_array_sub_a.change_basis(cb_op).map_to_asu()
      assert approx_equal(
        miller_array_sub_a.correlation(other=cb).coefficient(), 1)
    #
    coset_decomposition = sgtbx.cosets.left_decomposition_point_groups_only(
      g=lattice_group,
      h=subgroup_info.group())
    coset_decompositions.append(coset_decomposition)
    #
    for partition in coset_decomposition.partitions:
      s = partition[0]
      expected_match = s.is_unit_mx()
      cb_op = sgtbx.change_of_basis_op(s)
      cb = miller_array_sub_a.change_basis(cb_op).map_to_asu()
      is_match = abs(
        miller_array_sub_a.correlation(other=cb).coefficient() - 1) < 1.e-6
      assert is_match == expected_match
      for s in partition[1:]:
        cb_op = sgtbx.change_of_basis_op(s)
        cb = miller_array_sub_a.change_basis(cb_op).map_to_asu()
        assert approx_equal(cb.correlation(other=cb).coefficient(), 1)
    #
    blur_scale = 3
    blurred_data = miller_array_sub_a.data() * blur_scale
    # set to False to compare data with itself (will lead to failures below)
    if (True):
      blurred_data += flex.random_double(size=blurred_data.size())
    blurred_data /= blur_scale
    blurred_array = miller_array_sub_a.customized_copy(
      data=blurred_data)
    expected_correlation = "%.6f" % miller_array_sub_a.correlation(
      other=blurred_array).coefficient()
    if (verbose):
      print "blurred_array expected_correlation:", expected_correlation
    blurred_array_cb = blurred_array.change_basis(
        random.choice(list(subgroup_info.group())).as_xyz()).map_to_asu()
    self_ccs = {}
    for partition in coset_decomposition.partitions:
      cb_op = sgtbx.change_of_basis_op(partition[0])
      cb = miller_array_sub_a.change_basis(cb_op).map_to_asu()
      cc = blurred_array_cb.correlation(other=cb).coefficient()
      key = "%.6f" % cc
      self_ccs.setdefault(key, []).append(str(partition[0]))
    assert expected_correlation in self_ccs
    repetitions = [len(v) for v in self_ccs.values()]
    failure = (max(repetitions) > 1)
    if (failure or verbose):
      print subgroup_info
      print "self_ccs ops:", self_ccs.values()
    if (failure):
      for cc,ops in self_ccs.items():
        if (len(ops) > 1):
          print cc, ops
          for op in ops:
            ri = sgtbx.rt_mx(op).r().info()
            print "   ", ri.type(), ri.sense(), ri.ev()
      raise RuntimeError("max(repetitions) > 1")
  #
  for i in xrange(len(coset_decompositions)):
    for j in xrange(len(coset_decompositions)):
      exercise_double_coset_decomposition(
        crystal_symmetry_ri,
        lattice_group,
        miller_array_subs,
        miller_array_sub_as,
        miller_array_sub_bs,
        coset_decompositions,
        i, j,
        verbose)

def exercise_double_coset_decomposition(
      crystal_symmetry_ri,
      lattice_group,
      miller_array_subs,
      miller_array_sub_as,
      miller_array_sub_bs,
      coset_decompositions,
      i, j,
      verbose):
  group_a = miller_array_subs[i].space_group_info().group()
  group_b = miller_array_subs[j].space_group_info().group()
  miller_array_sub_a = miller_array_sub_as[i]
  miller_array_sub_b = miller_array_sub_bs[j]
  single_coset_ccs = {}
  for partition in coset_decompositions[i].partitions:
    cb_op = sgtbx.change_of_basis_op(partition[0])
    cb = miller_array_sub_a.change_basis(cb_op).map_to_asu()
    cc = cb.correlation(other=miller_array_sub_b).coefficient()
    key = "%.6f" % cc
    single_coset_ccs.setdefault(key, []).append(str(partition[0]))
  double_coset_ccs = {}
  for c in sgtbx.cosets.double_unique(lattice_group, group_a, group_b):
    cb_op = sgtbx.change_of_basis_op(c)
    cb = miller_array_sub_b.change_basis(cb_op).map_to_asu()
    cc = cb.correlation(other=miller_array_sub_a).coefficient()
    key = "%.6f" % cc
    double_coset_ccs.setdefault(key, []).append(str(c))
  double_coset_repetitions = [len(v) for v in double_coset_ccs.values()]
  failure = (max(double_coset_repetitions) > 1)
  if (failure or verbose):
    print [str(sgtbx.space_group_info(group=g)) for g in (group_a, group_b)]
    print "single_coset ops:", single_coset_ccs.values()
    print "double_coset ops:", double_coset_ccs.values()
  if (failure):
    for cc,ops in double_coset_ccs.items():
      if (len(ops) > 1):
        print cc, ops
        for op in ops:
          ri = sgtbx.rt_mx(op).r().info()
          print "   ", ri.type(), ri.sense(), ri.ev()
    raise RuntimeError("max(double_coset_repetitions) > 1")
  single_coset_ccs = single_coset_ccs.keys()
  double_coset_ccs = double_coset_ccs.keys()
  single_coset_ccs.sort()
  double_coset_ccs.sort()
  failure = (double_coset_ccs != single_coset_ccs)
  if (failure or verbose):
    print "single_coset_ccs:", single_coset_ccs
    print "double_coset_ccs:", double_coset_ccs
  if (failure):
    raise RuntimeError("double_coset_ccs != single_coset_ccs")

def run():
  verbose = "--verbose" in sys.argv[1:]
  quick = "--quick" in sys.argv[1:]
  harder = "--harder" in sys.argv[1:]
  hardest = "--hardest" in sys.argv[1:]
  for symbol in sgtbx.bravais_types.acentric:
    if (quick):
      if (symbol != "P 1 2 1"): continue
    elif (harder):
      if (symbol != "F 4 3 2"): continue
    elif (not hardest):
      if (symbol != "P 4 2 2"): continue
    space_group_info = sgtbx.space_group_info(symbol=symbol)
    for anomalous_flag in [False, True]:
      if (verbose):
        print symbol, "anomalous_flag =", anomalous_flag
      run_away_counter = 0
      while True:
        try:
          exercise(
            space_group_info=space_group_info,
            anomalous_flag=anomalous_flag,
            verbose=verbose)
        except RuntimeError, e:
          if (str(e) != "max(double_coset_repetitions) > 1"): raise
          print e, "(ignored since it may happen by chance)"
          run_away_counter += 1
          assert run_away_counter < 10
        else:
          break
  print "OK"

if (__name__ == "__main__"):
  run()
