from cctbx.development import random_structure
from cctbx.development.space_group_option_parser \
     import space_group_option_parser
from cctbx import miller
from cctbx import symmetry_search
from cctbx import sgtbx
from cctbx.array_family import flex
from scitbx import matrix as mat
import libtbx
import scitbx.random
import random

def exercise(space_group_info,
             fixed_random_seed=True,
             shifted_origin=None,
             elements=None,
             d_min=0.8,
             grid_resolution_factor=1/3.,
             verbose=False,
             **kwds):
  if elements is None:
    n_C = 5
    n_O = 1
    n_N = 1
    elements = ["C"]*n_C + ["O"]*n_O + ["N"]*n_N
  if verbose:
    print elements

  target_space_group_type = space_group_info.type()
  hall = sgtbx.space_group_symbols(
    target_space_group_type.lookup_symbol()).hall()
  print hall

  if fixed_random_seed:
    random.seed(1)
    flex.set_random_seed(1)

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
  if shifted_origin is None:
    shifted_origin = flex.random_double(3)
  shifted_origin = mat.col(shifted_origin)
  if verbose:
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
  if verbose:
    print sf_symm

  solution_hall, target_hall = [
    sgi.as_reference_setting().type().hall_symbol()
    for sgi in (sf_symm.space_group_info,
                target_structure.space_group_info()) ]
  assert solution_hall == target_hall, (solution_hall, target_hall)

  # Shift maximises goodness of symmetry?
  gos, solution_f = sf_symm.symmetrised_structure_factors()
  if space_group_info.type().hall_symbol() != ' P 1':
    assert gos.correlation > 0.99
    assert gos.gradient == (0, 0, 0)
    gos_away_from_max = sf_symm.symmetrised_structure_factors(
      delta=mat.col((0.1, 0.1, 0.1)))[0]
    assert gos_away_from_max.correlation < 0.9, gos_away_from_max.correlation

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
  ^           V solution_to_target_cb_op = (Q, q)
  ^           V
  ^------------

  The total transfrom from the target structure back to it reads
  (QP, q') with q' = (Q,q)(-s + P sigma) = (Q,q)delta_o
  with
  delta_o = sf_symm.cb_op_to_primitive(shifted_origin) - sf_symm.origin

  (Q, q') must leave the target structure space group invariant.
  Most of the time Q is the unit matrix and the test boils down to check
  whether delta is an allowed origin shift after changing to the input cell
  but it does not hurt to do the more general test all the time.
  """

  solution_to_target_cb_op = (
    target_structure.space_group_info()
    .change_of_basis_op_to_reference_setting().inverse()
    *
    sf_symm.space_group_info
    .change_of_basis_op_to_reference_setting())

  if verbose:
    print
    print "solution -> target: %s" % solution_to_target_cb_op.as_xyz()
  delta_o = (mat.col(sf_symm.cb_op_to_primitive(shifted_origin))
             - sf_symm.origin)
  delta = mat.col(solution_to_target_cb_op(delta_o))
  stabilising_cb_op = sgtbx.change_of_basis_op(sgtbx.rt_mx(
    (solution_to_target_cb_op*sf_symm.cb_op_to_primitive).c().r(),
    sgtbx.tr_vec((delta*72).as_int()*2, tr_den=144)))
    # guarding against rounding errors on some platforms (e.g. FC8)
    # meaning that (delta*144).as_int() would not work.
  target_sg = target_structure.space_group()
  assert target_sg == target_sg.change_basis(stabilising_cb_op)

def run():
  import sys
  libtbx.utils.show_times_at_exit()
  parser = space_group_option_parser()
  parser.option(None, '--skip_extra_tests',
                action='store_true',
                default=False)
  parser.option(None, '--fixed_random_seed',
                default=True)
  command_line = parser.process(sys.argv[1:])

  if not command_line.options.skip_extra_tests:
    # the solution-to-target change-of-basis computed as
    # reference-to-target x solution-to-reference
    # is not a mere translation
    exercise(
      sgtbx.space_group_info(
        'hall: -I 4 2c (1/2*x+1/2*y+1/12,-1/2*x+1/2*y-1/12,z-1/4)'),
      shifted_origin=(0, 0.4, 0.9),
      verbose=command_line.options.verbose)
    exercise(sgtbx.space_group_info('hall: C 2c 2 (x-y,x+y,z)'),
             shifted_origin=(0.4, 0.9, 0.6),
             verbose=command_line.options.verbose)

    # the centring translation search peaks (1/2, 0, 0) which is
    # a structure seminvariant: test whether this is rejected
    exercise(sgtbx.space_group_info('hall: -P 2a 2a'),
             shifted_origin=(0.1, 0.2, 0.6),
             verbose=command_line.options.verbose)

  # the traditional loop over a selection of space-groups
  # the last one, -I 4 2c (1/2*x+1/2*y+1/12,-1/2*x+1/2*y-1/12,z-1/4),
  # results in solution_to_target_cb_op being non-trivial.
  command_line.loop_over_space_groups(exercise,
                                      shifted_origin=(0.1, 0.2, 0.6))


if __name__ == '__main__':
  run()
