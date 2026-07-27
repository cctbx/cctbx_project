from __future__ import absolute_import, division, print_function

from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns_solving
from smtbx.refinement import cgls, least_squares
from smtbx.refinement.tests.tst_scipy_minimisers import (
  least_squares_for, observations_from, reference_structure, shaken_structure)
from libtbx.test_utils import approx_equal
import cmath
import math

# Weights which do not depend on Fc, so that the objective is a function of the
# parameters alone and the two paths can be held to agreeing exactly. What
# happens with the SHELX schemes is measured in tst_scipy_minimisers.
consistent_weighting = least_squares.mainstream_shelx_weighting(a=0, b=0)


def exercise_step_matches_cholesky():
  """ The CGLS step shall be the step the normal equations give.

  This is the whole claim of the module: that conjugate gradients on the
  design matrix solve the same system the Cholesky decomposition solves, not an
  approximation of it. The separable scale factor is what makes it worth
  checking -- the Gauss-Newton matrix is not J^T W J but D^T W D with
  D = k* J + yc (x) grad_k*, and getting that wrong would show up here.
  """
  observations = observations_from(reference_structure())
  for name, weighting in (('unit', least_squares.unit_weighting()),
                          ('sigma', least_squares.sigma_weighting()),
                          ('shelx a=0', consistent_weighting)):
    ls = least_squares_for(shaken_structure(), observations, weighting)
    ls.build_up()
    scale_factor = ls.scale_factor()
    ls.solve()
    reference = ls.step().deep_copy()

    other = least_squares_for(shaken_structure(), observations, weighting)
    other.build_up() # this is what makes the floating origin directions known
    problem = cgls.linearised_problem(other, scale_factor)
    step, iterations = cgls.solve(
      problem, preconditioner=problem.block_preconditioner(other),
      tolerance=1e-10)

    assert approx_equal(problem.scale_factor, scale_factor, eps=1e-12), name
    difference = (step - reference).norm()/reference.norm()
    assert difference < 1e-8, (name, difference)
    print("\t%-10s step agrees to %.1e in %d conjugate gradient iterations "
          "for %d parameters" % (name, difference, iterations, step.size()))


def exercise_block_preconditioner():
  """ The blocks shall be the diagonal blocks of the normal matrix.

  They are built without forming the normal matrix, by expanding D^T W D into
  terms already to hand, so they are worth checking against the matrix itself.
  """
  import numpy
  from libtbx import object_oriented_patterns as oop
  observations = observations_from(reference_structure())
  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)
  # The floating origin is projected out rather than added to the matrix, so
  # the blocks carry no trace of it and the matrix to compare them against must
  # not either. Only the Cholesky decomposition needs it, and nothing here
  # decomposes anything.
  ls.origin_fixing_restraint = oop.null()
  ls.build_up()
  n = ls.reparametrisation.n_independents
  a = ls.normal_matrix_packed_u().matrix_packed_u_as_symmetric()
  a = a.as_numpy_array().reshape(n, n)
  problem = cgls.linearised_problem(ls, ls.scale_factor())
  assert problem.deflation_basis is None

  worst = 0
  for indices, inverse in problem.block_preconditioner(ls):
    block = numpy.linalg.inv(inverse)
    expected = a[numpy.ix_(indices, indices)]
    worst = max(worst, numpy.abs(block - expected).max()
                       /numpy.abs(expected).max())
  print("\tblocks agree with the normal matrix to %.1e" % worst)
  assert worst < 1e-8, worst


def exercise_restraints():
  """ A restraints manager shall be handled, empty or not.

  A caller commonly passes one whether or not the structure has any restraints,
  so an empty manager is the ordinary case and not an edge case. The rows it
  contributes also arrive as a sparse matrix rather than a dense one, which is
  the other thing worth pinning here.
  """
  from cctbx import geometry_restraints
  from smtbx.refinement import restraints
  observations = observations_from(reference_structure())

  bond_proxies = geometry_restraints.shared_bond_simple_proxy()
  bond_proxies.append(
    geometry_restraints.bond_simple_proxy(
      i_seqs=(0, 1), distance_ideal=1.5, weight=10.))
  bond_proxies.append(
    geometry_restraints.bond_simple_proxy(
      i_seqs=(1, 2), distance_ideal=1.5, weight=10.))

  for label, proxies in (('empty', None), ('two bond restraints',
                                           bond_proxies)):
    xs = shaken_structure()
    if proxies is None:
      manager = restraints.manager()
    else:
      manager = restraints.manager(bond_proxies=proxies)
    ls = least_squares_for(xs, observations, consistent_weighting)
    ls.restraints_manager = manager
    ls.build_up()
    problem = cgls.linearised_problem(ls, ls.scale_factor())
    step, iterations = cgls.solve(
      problem, preconditioner=problem.block_preconditioner(ls),
      tolerance=1e-10)

    ls.solve()
    reference = ls.step()
    difference = (step - reference).norm()/reference.norm()
    assert difference < 1e-8, (label, difference)
    print("\t%-20s step agrees to %.1e" % (label, difference))


def exercise_f_mask():
  """ A solvent mask shall be part of the linearised problem.

  A mask changes Fc for every reflection, and omitting it does not fail: the
  conjugate gradients converge perfectly well on the wrong system and refine
  towards a structure which is not the one the objective is measured against.
  """
  reference = reference_structure()
  observations = observations_from(reference)
  # A stand-in for the solvent contribution: not physical, merely large enough
  # and varied enough that leaving it out cannot be mistaken for rounding.
  f_calc = reference.crystal_symmetry().build_miller_set(
    d_min=0.8, anomalous_flag=False).structure_factors_from_scatterers(
      reference, algorithm="direct").f_calc()
  f_mask = f_calc.customized_copy(data=flex.complex_double(
    [0.1*abs(f)*cmath.exp(1j*1.7*i) for i, f in enumerate(f_calc.data())]))

  # The system is compared rather than the step. A mask has fixed phases, so
  # unlike the model it does not move when the origin does, and it therefore
  # pins directions the origin fixing restraint assumes are free -- leaving the
  # two paths legitimately disagreeing about those directions by ~1e-4, the
  # restraint suppressing them and project() removing them outright. That is a
  # question about the origin, not about the mask, so the origin is taken out
  # of both and what the mask governs is checked exactly.
  import numpy
  from libtbx import object_oriented_patterns as oop

  def system_for(mask):
    ls = least_squares_for(shaken_structure(), observations,
                           consistent_weighting)
    ls.f_mask = mask
    ls.origin_fixing_restraint = oop.null()
    ls.build_up()
    return ls

  ls = system_for(f_mask)
  n = ls.reparametrisation.n_independents
  a = ls.normal_matrix_packed_u().matrix_packed_u_as_symmetric()
  a = a.as_numpy_array().reshape(n, n)
  b = ls.step_equations().right_hand_side()
  problem = cgls.linearised_problem(ls, ls.scale_factor())

  v = flex.double([math.sin(2.7*i) for i in range(n)])
  worst = 0
  for name, got, expected in (
      ('observables', problem.observables, ls.fc_sq.data()),
      ('right hand side', problem.right_hand_side(), b),
      ('A v', problem.normal_matrix_times(v),
       flex.double(a.dot(v.as_numpy_array()))),
      ('diag(A)', problem.normal_matrix_diagonal(),
       flex.double(numpy.diag(a).copy()))):
    difference = (got - expected).norm()/expected.norm()
    assert difference < 1e-12, (name, difference)
    worst = max(worst, difference)

  # and the mask has to have made a difference, or this proves nothing
  unmasked = system_for(None)
  moved = ((problem.right_hand_side()
            - cgls.linearised_problem(
                unmasked, unmasked.scale_factor()).right_hand_side()).norm()
           /b.norm())
  assert moved > 1e-2, moved
  print("\tsystem agrees to %.1e with the mask in, which moves it by %.1e"
        % (worst, moved))


def exercise_matrix_free_matches_stored():
  """ The two ways of applying A shall be the same way.

  The stored path keeps the design matrix; the matrix free one recomputes each
  row inside the product, so that nothing of size n_refl by n_par is ever
  allocated. They are the same algebra reached by different routes -- one in
  numpy over a stored array, one in C++ over the reflection loop -- and a
  system with two paths through it needs the rarely taken one held against the
  other or it rots. The stored path is itself pinned to the Cholesky step in
  exercise_step_matches_cholesky, so this chains onto that.
  """
  from cctbx import geometry_restraints
  from smtbx.refinement import restraints
  observations = observations_from(reference_structure())

  bond_proxies = geometry_restraints.shared_bond_simple_proxy()
  bond_proxies.append(
    geometry_restraints.bond_simple_proxy(
      i_seqs=(0, 1), distance_ideal=1.5, weight=10.))

  # weights which depend on Fc, restraints present, and a floating origin to
  # project: the awkward combination rather than the easy one
  for label, weighting, proxies in (
      ('unit', least_squares.unit_weighting(), None),
      ('shelx a=0.1', least_squares.mainstream_shelx_weighting(a=0.1, b=0.2),
       None),
      ('shelx + restraints',
       least_squares.mainstream_shelx_weighting(a=0.1, b=0.2), bond_proxies)):

    def problem_of(kind):
      ls = least_squares_for(shaken_structure(), observations, weighting)
      if proxies is not None:
        ls.restraints_manager = restraints.manager(bond_proxies=proxies)
      ls.build_up()
      return ls, kind(ls, ls.scale_factor())

    ls, stored = problem_of(cgls.linearised_problem)
    v = flex.double([math.sin(2.7*i) for i in range(stored.n_parameters)])
    step, iterations = cgls.solve(
      stored, preconditioner=stored.block_preconditioner(ls), tolerance=1e-10)

    for kind, sharp in ((cgls.matrix_free_problem, True),
                        (cgls.normal_matrix_problem, False)):
      other, alternative = problem_of(kind)
      assert approx_equal(alternative.scale_factor, stored.scale_factor,
                          eps=1e-12)
      # Both reach the same matrix by different routes -- one reproducing
      # D^T W D term by term in C++, one letting lstbx accumulate and finalise
      # it -- and both agree to rounding, so both are held to it.
      tolerance = 1e-10
      for name, got, expected in (
          ('right hand side', alternative.right_hand_side(),
           stored.right_hand_side()),
          ('A v', alternative.normal_matrix_times(stored.project(v)),
           stored.normal_matrix_times(stored.project(v)))):
        difference = (stored.project(got - expected).norm()
                      /stored.project(expected).norm())
        assert difference < tolerance, (label, kind, name, difference)

      other_step, other_iterations = cgls.solve(
        alternative, preconditioner=alternative.block_preconditioner(other),
        tolerance=1e-10)
      difference = (other_step - step).norm()/step.norm()
      assert difference < 1e-6, (label, kind, difference)
      print("\t%-18s %-22s step agrees to %.1e in %d iterations against %d"
            % (label, kind.__name__, difference, other_iterations, iterations))


def exercise_single_precision_design_matrix():
  """ A design matrix held in single precision shall give the same step.

  The store narrows; the products against it do not, accumulating in double
  through matrix_*_mixed. That is the whole claim, and it is what makes the
  halved matrix cost 2.6e-08 of relative error on a real structure rather than
  the 3.0e-06 a single precision accumulation would -- so what is checked here
  is not merely that it runs but that it lands far closer to double than a
  single precision accumulation could.

  Also checked: that the two builders really are different, since a wiring
  mistake which silently used the double one would pass every numerical test
  in this file.
  """
  observations = observations_from(reference_structure())
  ls = least_squares_for(shaken_structure(), observations,
                         least_squares.mainstream_shelx_weighting(a=0.1, b=0.2))
  ls.build_up()

  double = cgls.linearised_problem(ls, ls.scale_factor())
  single = cgls.linearised_problem(ls, ls.scale_factor(),
                                   single_precision=True)
  assert not double.built.is_mixed_precision()
  assert single.built.is_mixed_precision()
  assert not double.single_precision and single.single_precision

  # the products themselves, before any solve can average the error away
  v = flex.double([math.sin(2.7*i) for i in range(double.n_parameters)])
  u = flex.double([math.cos(1.3*i) for i in range(len(double.observables))])
  for name, got, expected in (
      ('J v', single._jacobian_times(v), double._jacobian_times(v)),
      ('J^T u', single._jacobian_transpose_times(u),
       double._jacobian_transpose_times(u)),
      ('right hand side', single.right_hand_side(), double.right_hand_side())):
    difference = (got - expected).norm()/expected.norm()
    # float has ~1.2e-07 of relative precision, so the storage rounding alone
    # puts a floor here; what this rules out is an accumulation in single,
    # which over these many terms would be orders worse
    assert difference < 1e-6, (name, difference)
    print("\t%-18s single precision store agrees to %.1e" % (name, difference))

  step_d, n_d = cgls.solve(
    double, preconditioner=double.block_preconditioner(ls), tolerance=1e-10)
  step_s, n_s = cgls.solve(
    single, preconditioner=single.block_preconditioner(ls), tolerance=1e-10)
  difference = (step_s - step_d).norm()/step_d.norm()
  assert difference < 1e-5, difference
  print("\tstep agrees to %.1e in %d conjugate gradient iterations against %d"
        % (difference, n_s, n_d))


def exercise_single_precision_refines():
  """ And a whole refinement with it shall reach the same place. """
  observations = observations_from(reference_structure())
  ends = {}
  for single in (False, True):
    ls = least_squares_for(shaken_structure(), observations,
                           least_squares.mainstream_shelx_weighting(a=0.1,
                                                                    b=0.2))
    steps = cgls.cgls_iterations(ls, n_max_iterations=5, mode='stored',
                                 single_precision_design_matrix=single)
    ends[single] = (ls.objective(), ls.scale_factor())
    # the reason is what a caller reports, so it has to say which precision
    # the budget was compared in
    assert ('single precision' in steps.mode_reason) == single, \
      steps.mode_reason
  difference = abs(ends[True][0] - ends[False][0])/abs(ends[False][0])
  assert difference < 1e-6, (ends, difference)
  print("\trefined objective %.8e single against %.8e double"
        % (ends[True][0], ends[False][0]))


def exercise_every_mode_refines():
  """ All three ways of holding the system shall refine to the same place.

  The objective is read after a build of the test's own rather than off the
  back of the run. The modes end their runs differently on purpose -- stored
  and matrix free close with an ordinary build, needing one for the covariance
  matrix, whereas the normal matrix mode already has one standing and would be
  paying a whole reflection loop purely to refresh a reported number. So what
  each leaves behind describes the structure at a different moment, and
  comparing those would measure when the last build happened rather than where
  the refinement arrived.
  """
  observations = observations_from(reference_structure())

  def converged(mode):
    ls = least_squares_for(shaken_structure(), observations,
                           consistent_weighting)
    cycles = cgls.cgls_iterations(ls, n_max_iterations=20, mode=mode)
    ls.build_up()
    return ls, cycles

  ls, stored = converged('stored')
  expected, expected_scale = ls.objective(), ls.scale_factor()

  for mode in ('matrix_free', 'normal_matrix'):
    other, cycles = converged(mode)
    assert approx_equal(other.objective(), expected, eps=1e-8), (
      mode, other.objective(), expected)
    assert approx_equal(other.scale_factor(), expected_scale, eps=1e-6), mode
    print("\t%-14s refined in %d cycles and %d conjugate gradient iterations "
          "to %.8e" % (mode, cycles.n_iterations, cycles.n_cg_iterations,
                       other.objective()))
  print("\t%-14s refined in %d cycles and %d conjugate gradient iterations "
        "to %.8e" % ('stored', stored.n_iterations, stored.n_cg_iterations,
                     expected))


def exercise_precomputed_mask_data():
  """ A caller which assembles the mask itself shall be honoured, not shadowed.

  A caller may build the MaskData once and keep it on the L.S. object as
  f_mask_data, the mask being fixed for the length of a run. A method of that
  same name on the base class would be shadowed by that attribute, and every
  path through build_up would fail, not merely this one.

  So the attribute is the cache and mask_data() is the accessor. Nothing here
  subclasses anything, which is exactly why this needs testing explicitly: the
  suite drives the base class, so a name a subclass happens to use is invisible
  to it.
  """
  from smtbx_refinement_least_squares_ext import MaskData
  reference = reference_structure()
  observations = observations_from(reference)
  f_calc = reference.crystal_symmetry().build_miller_set(
    d_min=0.8, anomalous_flag=False).structure_factors_from_scatterers(
      reference, algorithm="direct").f_calc()
  f_mask = f_calc.customized_copy(data=flex.complex_double(
    [0.1*abs(f)*cmath.exp(1j*1.7*i) for i, f in enumerate(f_calc.data())]))

  # the ordinary way, letting the L.S. object assemble it
  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)
  ls.f_mask = f_mask
  ls.build_up()
  expected = cgls.linearised_problem(ls, ls.scale_factor()).right_hand_side()

  # and the other way, handing over one already built. f_mask is deliberately
  # left unset, so anything which ignored the cache would silently compute an
  # unmasked system and the comparison would fail rather than pass vacuously.
  other = least_squares_for(shaken_structure(), observations,
                            consistent_weighting)
  other.f_mask_data = MaskData(
    other.observations, other.xray_structure.space_group(),
    other.observations.fo_sq.anomalous_flag(), f_mask.indices(), f_mask.data())
  assert other.f_mask is None
  other.build_up()
  got = cgls.linearised_problem(other, other.scale_factor()).right_hand_side()

  difference = (got - expected).norm()/expected.norm()
  assert difference < 1e-12, difference
  print("\ta mask handed over is used by build_up and by cgls alike (%.1e)"
        % difference)


def origin_fixed_structure(site_delta=0., u_delta=0.):
  """ The same atoms in P-1, where the symmetry fixes the origin.

  Wanted for the damping comparison. In P1 the origin floats, and the two paths
  treat that differently on purpose: lstbx restrains it and so damps it along
  with everything else, while cgls projects it out and has nothing there to
  damp. The restraint is a gauge fixing device rather than data, so damping it
  is meaningless either way -- but it does mean the two steps cannot be held
  against each other exactly while it is present. Where the origin is fixed the
  question does not arise.
  """
  from cctbx import crystal, xray
  from smtbx.refinement.tests import tst_scipy_minimisers as base
  cs = crystal.symmetry(unit_cell=base.unit_cell, space_group_symbol="P-1")
  xs = xray.structure(crystal.special_position_settings(cs))
  for i, site in enumerate(base.scatterer_sites):
    site = tuple(x + site_delta*math.sin(3.1*(3*i + j))
                 for j, x in enumerate(site))
    sc = xray.scatterer("C%i" % (i + 1), site,
                        u=base.u_isos[i] + u_delta*(-1)**i)
    sc.flags.set_use_u_iso(True).set_use_u_aniso(False)
    sc.flags.set_grad_site(True).set_grad_u_iso(True)
    xs.add_scatterer(sc)
  return xs


def exercise_damping():
  """ DAMP shall mean here what it means for Gauss-Newton.

  The first argument of ShelXL's DAMP multiplies the diagonal of the normal
  matrix by 1 + damp before the system is solved. lstbx does exactly that for
  the Cholesky path, in normal_eqns_solving.iterations.do_damping, and a caller
  feeds both paths from the same instruction -- so the same instruction has to
  produce the same step whichever solver is chosen, or it means two things.
  """
  observations = observations_from(origin_fixed_structure())
  shaken = lambda: origin_fixed_structure(site_delta=0.03, u_delta=0.004)

  for damping in (0.0007, 0.05):
    ls = least_squares_for(shaken(), observations, consistent_weighting)
    ls.build_up()
    assert not ls.origin_fixing_restraint.singular_directions
    # what lstbx does to the normal equations for this damping value
    a = ls.normal_matrix_packed_u()
    a.matrix_packed_u_diagonal_add_in_place(
      damping*a.matrix_packed_u_diagonal())
    ls.solve()
    expected = ls.step().deep_copy()

    for kind in (cgls.linearised_problem, cgls.normal_matrix_problem,
                 cgls.matrix_free_problem):
      other = least_squares_for(shaken(), observations, consistent_weighting)
      other.build_up()
      problem = kind(other, other.scale_factor())
      problem.set_damping(damping)
      step, iterations = cgls.solve(
        problem, preconditioner=problem.block_preconditioner(other),
        tolerance=1e-10)
      difference = (step - expected).norm()/expected.norm()
      assert difference < 1e-8, (damping, kind, difference)
    print("\tDAMP %-7g all three modes agree with the damped normal "
          "equations to %.1e" % (damping, difference))

  # and it has to be doing something, or the agreement proves nothing
  ls = least_squares_for(shaken(), observations, consistent_weighting)
  ls.build_up()
  problem = cgls.linearised_problem(ls, ls.scale_factor())
  plain, _ = cgls.solve(problem,
                        preconditioner=problem.block_preconditioner(ls),
                        tolerance=1e-10)
  problem.set_damping(0.05)
  damped, _ = cgls.solve(problem,
                         preconditioner=problem.block_preconditioner(ls),
                         tolerance=1e-10)
  moved = (damped - plain).norm()/plain.norm()
  assert moved > 1e-3, moved
  assert damped.norm() < plain.norm(), (damped.norm(), plain.norm())
  print("\tDAMP 0.05 shortens the step by %.1f%% and moves it by %.1e"
        % (100*(1 - damped.norm()/plain.norm()), moved))


def exercise_journal():
  """ Every mode shall leave a history of the run behind it.

  normal_eqns_solving decorates the L.S. object to record the objective, the
  gradient norm and the scale factor at each build, and a caller reads the last
  of those afterwards for the overall scale factor. A mode which builds by
  going straight to the undecorated object records nothing, and since this one
  does the only builds there are, that would be the entire history -- an empty
  list, and an IndexError in the caller rather than anywhere near here.
  """
  observations = observations_from(reference_structure())
  for mode in ('stored', 'normal_matrix', 'matrix_free'):
    ls = least_squares_for(shaken_structure(), observations,
                           consistent_weighting)
    cycles = cgls.cgls_iterations(ls, n_max_iterations=20, mode=mode)
    history = cycles.non_linear_ls.journal
    assert history.scale_factor_history.size() > 0, mode
    assert history.objective_history.size() > 0, mode
    assert (history.scale_factor_history.size()
            == history.objective_history.size()), mode
    # and the last entry describes the model the run ended on, not the scale
    # factor an earlier cycle was weighted with
    assert approx_equal(history.scale_factor_history[-1],
                        ls.optimal_scale_factor(), eps=1e-8), mode
    print("\t%-14s journalled %d builds, closing scale factor %.6f"
          % (mode, history.scale_factor_history.size(),
             history.scale_factor_history[-1]))


def exercise_mode_selection():
  """ The mode shall be chosen on what fits, largest matrix first. """
  observations = observations_from(reference_structure())
  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)
  ls.build_up()
  chooser = cgls.cgls_iterations.__new__(cgls.cgls_iterations)
  chooser.mode = None

  chooser.max_design_matrix_memory = 1e9
  chooser.max_normal_matrix_memory = 1e9
  assert chooser.mode_for(ls) == 'stored'
  # the design matrix no longer fits, but the normal matrix does
  chooser.max_design_matrix_memory = 0
  assert chooser.mode_for(ls) == 'normal_matrix'
  # nor does that
  chooser.max_normal_matrix_memory = 0
  assert chooser.mode_for(ls) == 'matrix_free'
  print("\tmode falls back stored -> normal_matrix -> matrix_free")

  # And the default, which is to ask the machine rather than to hold a number.
  # This is the path that matters: a fixed 2048 MB default once excluded a
  # 2555 MB design matrix on a host with 128 GB, and sent a structure to the
  # normal matrix at three times the cost a cycle. Nothing above would have
  # caught it, every case there pinning the budget by hand.
  chooser.max_design_matrix_memory = None
  chooser.max_normal_matrix_memory = 1e9
  budget = chooser.design_matrix_budget()
  assert budget > 0, budget
  assert chooser.mode_for(ls) == 'stored'
  # a share of the machine, not the whole of it
  from libtbx.utils import guess_total_memory
  total = guess_total_memory()/1048576.
  assert budget < total, (budget, total)
  print("\tby default the budget is %.0f MB of %.0f MB on this machine"
        % (budget, total))


def exercise_refinement():
  """ A CGLS refinement shall arrive where Gauss-Newton arrives. """
  observations = observations_from(reference_structure())

  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)
  normal_eqns_solving.naive_iterations(
    ls, gradient_threshold=1e-12, step_threshold=1e-7)
  best, best_scale = ls.objective(), ls.scale_factor()

  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)
  cycles = cgls.cgls_iterations(ls, n_max_iterations=20)
  assert cycles.n_iterations < 20, "did not converge"
  assert approx_equal(ls.objective(), best, eps=1e-10), (ls.objective(), best)
  assert approx_equal(ls.scale_factor(), best_scale, eps=1e-8)
  print("\trefined in %d cycles and %d conjugate gradient iterations, "
        "objective %.8e against Gauss-Newton's %.8e"
        % (cycles.n_iterations, cycles.n_cg_iterations, ls.objective(), best))


def exercise_standard_uncertainties():
  """ The covariance matrix shall be available once a CGLS run is over.

  CGLS cannot produce one, there being no normal matrix to invert, so the run
  has to end with normal equations standing. Anyone refining a structure wants
  s.u. afterwards, so this is not optional.

  Every mode is exercised because they arrive at that differently: stored and
  matrix free close with a build of their own, while the normal matrix mode
  relies on the last cycle's still being there and skips the closing build
  entirely. That saved a whole reflection loop, and it is exactly the sort of
  saving that would go unnoticed until someone asked for an s.u. and got a
  singular matrix.
  """
  observations = observations_from(reference_structure())
  for mode in ('stored', 'normal_matrix', 'matrix_free'):
    ls = least_squares_for(shaken_structure(), observations,
                           consistent_weighting)
    cgls.cgls_iterations(ls, n_max_iterations=20, mode=mode)
    assert ls.overridden_scale_factor is None, mode
    variances = ls.covariance_matrix_and_annotations().diagonal()
    assert variances.size() == ls.reparametrisation.n_independents, mode
    for i, variance in enumerate(variances):
      assert variance > 0 and not math.isnan(variance), (mode, i)
    print("\t%-14s s.u. span %.2e to %.2e"
          % (mode, math.sqrt(flex.min(variances)),
             math.sqrt(flex.max(variances))))

  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)
  cgls.cgls_iterations(ls, n_max_iterations=20)
  # and the override the run used must not be left behind
  assert ls.overridden_scale_factor is None
  variances = ls.covariance_matrix_and_annotations().diagonal()
  assert variances.size() == ls.reparametrisation.n_independents
  for i, variance in enumerate(variances):
    assert variance > 0 and not math.isnan(variance), i
  print("\tstandard uncertainties span %.2e to %.2e"
        % (math.sqrt(flex.min(variances)), math.sqrt(flex.max(variances))))


def exercise_without_standard_uncertainties():
  """ Declining the s.u. shall cost the s.u. and nothing else.

  The closing build is a full Gauss-Newton one and costs about three cycles of
  this method, so a refinement whose result is another refinement can decline
  it. What it must not cost is the report: report_state describes the model
  each cycle *entered* with, running before the step it reports on is applied,
  so a run which simply stopped would leave the objective and the structure
  factors one shift out of date and ls.objective() still answering from the
  opening build. The run closes on a gradient-free pass instead, and this holds
  it to the same answer as the full one.
  """
  observations = observations_from(reference_structure())
  refined = {}
  for wanted in (True, False):
    ls = least_squares_for(shaken_structure(), observations,
                           consistent_weighting)

    class iterations(cgls.cgls_iterations):
      compute_standard_uncertainties = wanted

    iterations(ls, n_max_iterations=20, mode='stored')
    assert ls.overridden_scale_factor is None, wanted
    refined[wanted] = (ls.objective(), ls.objective_data_only,
                       ls.optimal_scale_factor())
  for i, name in enumerate(('objective', 'objective_data_only',
                            'scale factor')):
    assert approx_equal(refined[True][i], refined[False][i], eps=1e-12), name
  print("\tdeclining the s.u. leaves the objective at %.8e either way"
        % refined[False][0])

  # and the s.u. really are the only thing given up
  ls = least_squares_for(shaken_structure(), observations,
                         consistent_weighting)

  class no_su(cgls.cgls_iterations):
    compute_standard_uncertainties = False

  no_su(ls, n_max_iterations=20, mode='stored')
  try:
    ls.covariance_matrix_and_annotations()
  except Exception:
    pass          # no normal matrix to invert, which is the whole point
  else:
    # a matrix left over from an earlier build would be worse than none: it
    # would describe a model the run has since moved away from
    raise AssertionError(
      "a covariance matrix was available after declining to build one")


def exercise_no_flag_shadows_a_method():
  """ A flag adopted onto the instance shall not shadow a method of the class.

  Olex2 mixes cgls_iterations into a class that carries a standard_uncertainties
  *method* -- scale_shifts calls it, per cycle, to damp the step by shift/su.
  Every optional flag on cgls_iterations reaches the instance through
  adopt_optional_init_args, which setattrs it; so a flag sharing a method's name
  lands on the instance as a bool and the next call to that method finds the
  bool instead. That is not a hypothetical: a flag named `standard_uncertainties`
  did exactly this and broke every CGLS-J refinement run from Olex2, invisibly
  to a test suite that never mixes such a method in. So mix one in here.

  This is a naming contract, checked structurally: no class attribute of
  cgls_iterations whose value is a plain bool may share a name with a method it
  could be mixed in alongside. The concrete reproduction below is the reason;
  the loop is the guard against the next one.
  """
  class carries_a_method(cgls.cgls_iterations):
    # stand-in for iterations_with_shift_analysis.standard_uncertainties
    def standard_uncertainties(self):
      return None

  observations = observations_from(reference_structure())
  ls = least_squares_for(shaken_structure(), observations,
                         least_squares.mainstream_shelx_weighting(a=0.1, b=0.2))
  # pass every bool flag as a keyword, exactly as refinement.py does, so any
  # one of them shadowing the method would surface here
  steps = carries_a_method(ls, n_max_iterations=2, mode='stored',
                           compute_standard_uncertainties=False,
                           single_precision_design_matrix=True)
  assert callable(steps.standard_uncertainties), \
    "a flag shadowed the standard_uncertainties method"
  assert steps.standard_uncertainties() is None

  # and the general rule, so the next flag added cannot reintroduce this
  import types
  for name in dir(cgls.cgls_iterations):
    if name.startswith('__'):
      continue
    value = getattr(cgls.cgls_iterations, name)
    if isinstance(value, bool):
      # nothing this object could be mixed in alongside may define `name` as a
      # method; the two olex2 iteration classes are the live example
      assert not isinstance(
        getattr(carries_a_method, name, None), types.FunctionType), \
        ("flag %r shares its name with a method it is mixed in beside" % name)
  print("\tno bool flag shadows a method it is mixed in beside")


def run():
  exercise_step_matches_cholesky()
  exercise_block_preconditioner()
  exercise_restraints()
  exercise_f_mask()
  exercise_precomputed_mask_data()
  exercise_matrix_free_matches_stored()
  exercise_single_precision_design_matrix()
  exercise_single_precision_refines()
  exercise_no_flag_shadows_a_method()
  exercise_every_mode_refines()
  exercise_damping()
  exercise_journal()
  exercise_mode_selection()
  exercise_refinement()
  exercise_standard_uncertainties()
  exercise_without_standard_uncertainties()
  print('OK')


if __name__ == '__main__':
  run()
