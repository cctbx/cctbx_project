from cctbx.development import debug_utils, random_structure
from cctbx import miller
from cctbx import symmetry_search
from cctbx import sgtbx, maptbx
from cctbx import euclidean_model_matching as emma
from cctbx.array_family import flex
from scitbx import matrix as mat
import libtbx
from libtbx.test_utils import approx_equal
import scitbx.random
import random
import math

def exercise_one(flags, space_group_info,
                 shifted_origin=(0.1, 0.2, 0.6),
                 elements=None,
                 d_min=0.8,
                 grid_resolution_factor=1/3):
  if elements is None:
    n_C = 5
    n_O = 1
    n_N = 1
    elements = ["C"]*n_C + ["O"]*n_O + ["N"]*n_N
  if flags.Verbose:
    print elements

  target_space_group_type = space_group_info.type()
  hall = sgtbx.space_group_symbols(
    target_space_group_type.lookup_symbol()).hall()
  print hall

  # Generate a random structure in real space,
  # compute its structure factors,
  # that we will try to recover the symmetry of
  target_structure = random_structure.xray_structure(
    space_group_info=space_group_info,
    elements=elements,
    use_u_iso=True,
    random_u_iso=True,
    random_u_iso_scale=0.04,
    use_u_aniso=False,
  )
  if flags.non_random_shift:
    random.seed(1)
    flex.set_random_seed(1)
  else:
    shifted_origin = flex.random_double(3)
  shifted_origin = mat.col(shifted_origin)
  if flags.Verbose:
    print "new origin = (%.3f, %.3f, %.3f)" % shifted_origin.elems
    print
  target_structure_in_p1 = target_structure\
                         .expand_to_p1().apply_shift(shifted_origin)
  target_f_in_p1 = miller.build_set(
    crystal_symmetry=target_structure_in_p1,
    anomalous_flag=False,
    d_min=d_min
    ).structure_factors_from_scatterers(
      xray_structure=target_structure_in_p1,
      algorithm="direct").f_calc()

  # Recover space group?
  sf_symm = symmetry_search.structure_factor_symmetry(target_f_in_p1)
  assert (
    sf_symm.space_group_info_in_input_cell.group()\
    .conventional_centring_type_symbol()
    == target_structure.space_group().conventional_centring_type_symbol())
  if flags.Verbose:
    print sf_symm

  solution_hall, target_hall = [
    sgi.as_reference_setting().type().hall_symbol()
    for sgi in (sf_symm.space_group_info_in_input_cell,
                target_structure.space_group_info()) ]
  assert solution_hall == target_hall, (solution_hall, target_hall)

  # Shift maximises goodness of symmetry?
  gos, solution_f = sf_symm.symmetrised_structure_factors_in_input_cell()
  if space_group_info.type().hall_symbol() != ' P 1':
    assert gos.correlation > 0.99
    assert gos.gradient == (0, 0, 0)
    assert sf_symm.symmetrised_structure_factors_in_input_cell(
      delta=mat.col((0.1, 0.1, 0.1)))[0].correlation < 0.9

  # Recovered origin
  """The sequence of transform read:
  ----->target structure
  ^           V
  ^           V  (1, shifted_origin=sigma)
  ^           V
  ^      shifted target
  ^           V
  ^           V sf_symm.cb_op_to_primitive = (P, 0)
  ^           V
  ^      shifted target in primitive cell = shifted solution in primitive cell
  ^           V
  ^           V (1, -sf_symm.origin=-s)
  ^           V
  ^      solution in primitive cell
  ^           V
  ^           V (sf_symm.cb_op_to_primitive.inverse()) = (P^{-1}, 0)
  ^           V
  ^      solution in input cell
  ^           V
  ^           V solution_to_target_cb_op = (Q, q)
  ^           V
  ^------------

  The total transfrom from the target structure back to it reads
  (Q, q') with q' = (Q,q)(-P^{-1}s + sigma) = (Q,q)delta_o
  with
  delta_o = shifted_origin - sf_symm.origin_in_input_cell

  (Q, q') must leave the target structure space group invariant.
  Most of the time Q is the unit matrix and the test boils down to check
  whether delta is an allowed origin shift but it does not hurt to
  do the more general test all the time.
  """

  solution_to_target_cb_op = (
    target_structure.space_group_info()
    .change_of_basis_op_to_reference_setting().inverse()
    *
    sf_symm.space_group_info_in_input_cell
    .change_of_basis_op_to_reference_setting())

  if flags.Verbose:
    print
    print "solution -> target: %s" % solution_to_target_cb_op.as_xyz()
  delta_o = shifted_origin - sf_symm.origin_in_input_cell
  delta = mat.col(solution_to_target_cb_op(delta_o))
  stabilising_cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(
    solution_to_target_cb_op.c().r(),
    sgtbx.tr_vec((delta*144).as_int(), tr_den=144)))
  target_sg = target_structure.space_group()
  assert target_sg == target_sg.change_basis(stabilising_cb_op)

def exercise(argv):
  debug_utils.parse_options_loop_space_groups(
    argv,
    keywords=('non_random_shift',),
    call_back=exercise_one,
    symbols_to_stderr=False,
    d_min=0.8,
    grid_resolution_factor=1/3,
  )

def run():
  import sys, os
  if 0 and 'WINGDB_ACTIVE' in os.environ:
    sys.argv[1:] = [
      #'--non_random_shift',
      '--Verbose',
      #'hall: -I 4 2c (1/2*x+1/2*y+1/12,-1/2*x+1/2*y-1/12,z-1/4)',
      'hall: C 2c 2 (x+y, x-y, -z)'
      ]
    for i in xrange(100):
      exercise(sys.argv[1:])
    return

  if not sys.argv[1:]:
    flags = libtbx.group_args(Verbose='--Verbose' in sys.argv[1:],
                              non_random_shift=True)

    # the solution-to-target change-of-basis computed as
    # reference-to-target x solution-to-reference
    # is not a mere translation
    exercise_one(
      flags,
      sgtbx.space_group_info(
        'hall: -I 4 2c (1/2*x+1/2*y+1/12,-1/2*x+1/2*y-1/12,z-1/4)'),
      shifted_origin=(0, 0.4, 0.9))
    exercise_one(flags, sgtbx.space_group_info('hall: C 2c 2 (x-y,x+y,z)'),
                 shifted_origin=(0.4, 0.9, 0.6))

    # the centring translation search peaks (1/2, 0, 0) which is
    # a structure seminvariant: test whether this is rejected
    exercise_one(flags, sgtbx.space_group_info('hall: -P 2a 2a'))

  # the traditional loop over a selection of space-groups
  # the last one, -I 4 2c (1/2*x+1/2*y+1/12,-1/2*x+1/2*y-1/12,z-1/4),
  # results in solution_to_target_cb_op being non-trivial.
  exercise(sys.argv[1:])


if __name__ == '__main__':
  run()
