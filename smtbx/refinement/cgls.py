""" Conjugate gradient least squares, after Hendrickson & Konnert.

The refinement engines in smtbx form the normal matrix and solve it by Cholesky
decomposition. That costs O(n_refl n_par^2) to accumulate and O(n_par^3) to
factorise, both of which come to dominate as the number of parameters grows.
CGLS instead iterates conjugate gradients on the same system using only
matrix-vector products, so the normal matrix is never formed and never
factorised. It is the approach ShelXL takes with its CGLS instruction.

What is solved is exactly the system the Cholesky path solves, not an
approximation of it. That is worth spelling out, because the separable scale
factor makes the Gauss-Newton matrix something other than the obvious J^T W J.
Writing yc for the computed observables, yo for the observed, and k* for the
scale factor which the objective has optimised away,

    k*      = (w yo . yc) / (w yc . yc)
    r       = yo - k* yc
    grad_k* = (J^T W r - k* J^T W yc) / (w yc . yc)
    D       = k* J + yc (x) grad_k*

and then A = D^T W D and b = D^T W r, which is what
non_linear_ls_with_separable_scale_factor::finalise builds. So a product with
D costs one product with J plus a rank-one correction, and the whole solve
needs J only through J v and J^T u.

Restraints and the floating origin restraint enter as further rows of D, with
their own residuals and weights, exactly as they enter the normal equations.

The one thing CGLS cannot give is a covariance matrix, since that comes from
inverting the normal matrix. A refinement which needs standard uncertainties --
which in practice means any refinement whose results are to be published -- has
to end with one ordinary build and solve. ShelXL has the same limitation and
the same remedy.
"""
from __future__ import absolute_import, division, print_function

from cctbx import xray
from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns_solving
from smtbx_refinement_least_squares_ext import (
  build_design_matrix, build_design_matrix_single)
import sys


class linearised_problem(object):
  """ The Gauss-Newton system at one point, held as rows rather than a matrix.

  Built from a crystallographic_ls whose reparametrisation is up to date. The
  scale factor used to compute the weights has to be supplied, since the SHELX
  weighting schemes depend on it; passing the value from the previous cycle is
  what build_up does.

  The solvent contribution is taken from the L.S. object, exactly as build_up
  takes it. It has to be: without it the structure factors are not the ones the
  refinement is minimising against, and the whole cycle -- step, objective and
  reported residuals alike -- describes a different structure.
  """

  # How many bytes one element of a stored design matrix occupies, in each
  # precision, for the size estimate mode_for makes before building anything.
  bytes_per_element = {False: 8, True: 4}

  def __init__(self, ls, scale_factor, f_mask_data=None,
               single_precision=False):
    from smtbx_refinement_least_squares_ext import (
      separable_scale_factor_summary)
    reparametrisation = ls.reparametrisation
    fc_correction = reparametrisation.fc_correction
    if fc_correction is None:
      fc_correction = xray.dummy_fc_correction()
    if f_mask_data is None:
      f_mask_data = ls.mask_data()
    jacobian_transpose = \
      reparametrisation.jacobian_transpose_matching_grad_fc()
    self.n_parameters = jacobian_transpose.n_rows

    summary = self._summary_accumulator(ls)

    # One pass builds J and, in the same pass, the summary over it: the scale
    # factor, the right hand side and the preconditioner blocks. Recovering
    # those from J afterwards would be a second pass over the whole matrix.
    #
    # The threading flags are passed exactly as build_up passes them, since they
    # default off and leaving them out builds the matrix single-threaded.
    # use_openmp is not among them: that path forms a normal matrix per thread
    # and so cannot feed the accumulator. See least_squares.h.
    #
    # single_precision narrows the store, not the arithmetic: the products still
    # accumulate in double, so the matrix is half the size for a negligible
    # error. See matrix_vector_mixed in scitbx/matrix/matrix_vector_operations.h.
    builder = build_design_matrix_single if single_precision \
              else build_design_matrix
    self.single_precision = single_precision
    built = builder(
      summary, ls.observations, f_mask_data, ls.weighting_scheme, scale_factor,
      ls.one_h_linearisation, jacobian_transpose, fc_correction, False,
      getattr(ls, 'may_parallelise', False), False,
      getattr(ls, 'max_memory', 300))

    # The builder keeps the design matrix and does the products against it, so
    # nothing here holds a copy of it.
    self.built = built
    self.observables = flex.double(built.observables())
    self.weights = flex.double(built.weights())
    self.f_calc = built.f_calc()

    # The builder stored the rows already multiplied by sqrt(w), which is what
    # lets the matrix-vector products be plain BLAS on them. Only the inverse
    # is wanted here, to undo that where a product is taken against unweighted
    # quantities; a zero weight has no inverse and cannot contribute, so it is
    # given one that leaves the row alone.
    self.sqrt_weights = flex.sqrt(self.weights)
    safe = self.sqrt_weights.deep_copy()
    safe.set_selected(safe <= 0, 1.)
    self.inverse_sqrt_weights = 1./safe

    self.sum_w_yo_sq = summary.sum_w_yo_sq()
    self.weighted_observables_sq = summary.sum_w_yc_sq()
    self.scale_factor = summary.scale_factor()
    self.grad_scale_factor = summary.grad_scale_factor()
    self.residuals = (ls.observations.fo_sq.data()
                      - self.scale_factor*self.observables)
    self._right_hand_side = summary.right_hand_side()
    self._blocks = summary.blocks().as_numpy_array()

    self.extra_design_matrix = None
    self.extra_residuals = None
    self.extra_weights = None
    self.deflation_basis = None
    self._add_restraints(ls)
    self._add_floating_origin(ls)

  def objective(self):
    """ The value the refinement minimises.

    Normalised by sum w yo^2 as lstbx does, but note that only the structure
    factor equations are: finalise() normalises those, and the restraints and
    the floating origin restraint are added to the normal equations afterwards
    and so enter unscaled. Everything here follows the same convention, so that
    what comes out is the very system the Cholesky path solves.
    """
    residual = flex.sum(self.weights*flex.pow2(self.residuals))
    residual /= self.sum_w_yo_sq
    if self.extra_residuals is not None:
      residual += flex.sum(
        self.extra_weights*flex.pow2(self.extra_residuals))
    return residual/2

  def _summary_accumulator(self, ls):
    """ What the reflection pass accumulates besides whatever it stores.

    The blocks have to be settled before the pass which fills them, the
    accumulator gathering each block's own outer product as it goes, so this
    has to happen first and not on demand afterwards.
    """
    from smtbx_refinement_least_squares_ext import (
      separable_scale_factor_summary)
    self.blocks = self.parameter_blocks(ls)
    self.block_parameters = flex.int([int(i) for b in self.blocks for i in b])
    self.block_sizes = flex.int([len(b) for b in self.blocks])
    return separable_scale_factor_summary(
      self.n_parameters, self.block_parameters, self.block_sizes)

  def _jacobian_times(self, v):
    """ J v, undoing the sqrt(w) the stored matrix carries. """
    product = self.built.times(v)
    product *= self.inverse_sqrt_weights
    return product

  def _jacobian_transpose_times(self, u):
    """ J^T u, likewise. """
    return self.built.transpose_times(u*self.inverse_sqrt_weights)

  def _add_rows(self, design_matrix, residuals, weights):
    """ Append equations which are not structure factor equations.

    There are only ever a handful of these -- the restraints and the two or
    three floating origin directions -- so they are kept as one small dense
    block rather than being folded into the design matrix of the reflections.
    """
    import numpy
    # the restraints arrive as a sparse matrix and the origin directions as a
    # dense one, and neither is a numpy array to begin with
    if hasattr(design_matrix, 'as_dense_matrix'):
      design_matrix = design_matrix.as_dense_matrix()
    if hasattr(design_matrix, 'as_numpy_array'):
      design_matrix = design_matrix.as_numpy_array()
    rows = numpy.asarray(design_matrix).reshape(-1, self.n_parameters)
    if rows.shape[0] == 0: return
    residuals = flex.double(residuals)
    weights = flex.double(weights)
    if self.extra_design_matrix is None:
      self.extra_design_matrix = rows
      self.extra_residuals = residuals
      self.extra_weights = weights
    else:
      self.extra_design_matrix = numpy.vstack(
        (self.extra_design_matrix, rows))
      self.extra_residuals.extend(residuals)
      self.extra_weights.extend(weights)

  def _add_restraints(self, ls):
    # A refinement usually has a restraints manager whether or not there is
    # anything in it, so an empty one is the ordinary case, not an edge case.
    if ls.restraints_manager is None: return
    linearised = ls.restraints_manager.build_linearised_eqns(
      ls.xray_structure, ls.reparametrisation.parameter_map())
    if linearised.n_restraints() == 0: return
    jacobian = ls.reparametrisation.jacobian_transpose_matching(
      ls.reparametrisation.mapping_to_grad_fc_all).transpose()
    # the normalisation which puts the restraints on the scale of the average
    # residual, as build_up does
    factor = ls.restraints_normalisation_factor
    if factor is None:
      factor = 1
    # non_linear_ls::add_equations negates the right hand side by default, so
    # what build_up puts into the normal equations for a restraint is
    # -w delta G. The residual is stored negated here to say the same thing,
    # which leaves the objective alone, that going as delta squared.
    self._add_rows(linearised.design_matrix*jacobian,
                   -linearised.deltas,
                   linearised.weights*factor)

  def _add_floating_origin(self, ls):
    """ The directions along which the structure is free to drift.

    These are equations with a residual of zero: they add nothing to the right
    hand side and exist only to stop the otherwise singular normal matrix from
    defeating the Cholesky decomposition. They are not added to the system here,
    because for conjugate gradients they would be ruinous. Their weights are
    enormous beside those of the data -- two orders of magnitude or more -- and
    conjugate gradients converge at a rate set by the condition number, which
    such equations raise by as much again.

    What the restraint is trying to express is that the step has no component
    along those directions, and for an iterative solver that is said directly,
    by projecting them out of the search space. The directions are kept here as
    an orthonormal basis for project() to use.

    They are only known once the restraint has been applied to a set of normal
    equations at least once, which the initial build does.
    """
    import numpy
    directions = ls.origin_fixing_restraint.singular_directions
    if not directions: return
    jacobian_transpose = \
      ls.reparametrisation.jacobian_transpose_matching_grad_fc()
    rows = numpy.array(
      [(jacobian_transpose*flex.double(d)).as_numpy_array()
       for d in directions])
    basis, _ = numpy.linalg.qr(rows.T)
    self.deflation_basis = basis.T

  def project(self, v):
    """ Strip out any component along the directions the structure may drift.
    """
    if self.deflation_basis is None: return v
    q = self.deflation_basis
    a = v.as_numpy_array()
    return flex.double(a - q.T.dot(q.dot(a)))

  # The first argument of ShelXL's DAMP: the diagonal of the normal matrix is
  # multiplied by 1 + damp before the system is solved, which shortens the step
  # along the directions the data determine least well. It has to mean here what
  # normal_eqns_solving.iterations.do_damping makes it mean for the Cholesky
  # path, or the same instruction refines a structure two ways.
  damping = 0.
  damping_diagonal = None

  def set_damping(self, value):
    """ Damp the diagonal of A by 1 + value, as ShelXL's DAMP asks.

    The diagonal here is that of the data and the restraints. Where the origin
    floats it differs from the one the Cholesky path damps, which carries the
    origin fixing restraint as well -- that being a gauge fixing device rather
    than data, projected out here rather than restrained. Damping it is
    meaningless on either reading, but it does mean the two steps agree exactly
    only where the symmetry fixes the origin, and otherwise differ by a little,
    in proportion to the damping.
    """
    self.damping = value
    if value:
      self.damping_diagonal = self.normal_matrix_diagonal()
    else:
      self.damping_diagonal = None

  def normal_matrix_times(self, v):
    result = self.undamped_normal_matrix_times(v)
    if self.damping:
      result += self.damping*self.damping_diagonal*v
    return result

  def block_preconditioner(self, ls, max_block_size=24):
    """ The diagonal blocks of A, damped as A itself is, and inverted. """
    import numpy
    inverses = []
    for indices, block in self.raw_blocks(ls, max_block_size):
      if self.damping:
        diagonal = self.damping_diagonal
        for j, i in enumerate(indices):
          block[j, j] += self.damping*diagonal[int(i)]
      try:
        inverse = numpy.linalg.inv(block)
      except numpy.linalg.LinAlgError:
        # a block the data does not determine; fall back on its diagonal
        diagonal = numpy.diag(block).copy()
        diagonal[diagonal <= 0] = 1.
        inverse = numpy.diag(1./diagonal)
      inverses.append((indices, inverse))
    return inverses

  def undamped_normal_matrix_times(self, v):
    """ A v, with A = D^T W D, without ever forming A. """
    # data equations: D v = k* (J v) + (grad_k* . v) yc
    jv = self._jacobian_times(v)
    jv *= self.scale_factor
    jv += self.observables*flex.sum(self.grad_scale_factor*v)
    jv *= self.weights
    # D^T (W D v) = k* J^T u + (yc . u) grad_k*
    result = self._jacobian_transpose_times(jv)
    result *= self.scale_factor
    result += self.grad_scale_factor*flex.sum(self.observables*jv)
    result /= self.sum_w_yo_sq
    if self.extra_design_matrix is not None:
      rows = self.extra_design_matrix
      u = rows.dot(v.as_numpy_array())*self.extra_weights.as_numpy_array()
      result += flex.double(rows.T.dot(u))
    return result

  def right_hand_side(self):
    """ b = D^T W r, as the pass which built the system already worked it out.
    """
    result = self._right_hand_side.deep_copy()
    if self.extra_design_matrix is not None:
      rows = self.extra_design_matrix
      u = (self.extra_weights*self.extra_residuals).as_numpy_array()
      result += flex.double(rows.T.dot(u))
    return result

  def parameter_blocks(self, ls, max_block_size=24):
    """ The independent parameters grouped by the scatterer they belong to.

    A crystallographic normal matrix is strongly correlated within one atom --
    its coordinates against its ADPs -- and much more weakly between atoms, so
    these blocks capture most of what makes the matrix ill conditioned.

    A parameter which drives more than one scatterer, as constrained groups and
    riding hydrogens do, is put with the first of them: the grouping only has
    to be good enough to precondition on, not exact.
    """
    import numpy
    reparametrisation = ls.reparametrisation
    jacobian = \
      reparametrisation.jacobian_transpose_matching_grad_fc().transpose()
    annotations = reparametrisation.component_annotations
    groups, order = {}, []
    for i in range(self.n_parameters):
      label = None
      for component in jacobian.col(i):
        label = annotations[component[0]].split('.')[0]
        break
      if label not in groups:
        groups[label] = []
        order.append(label)
      groups[label].append(i)
    blocks = []
    for label in order:
      indices = groups[label]
      # a very large group would cost more to invert than it is worth
      for start in range(0, len(indices), max_block_size):
        blocks.append(numpy.array(indices[start:start + max_block_size]))
    return blocks

  def raw_blocks(self, ls, max_block_size=24):
    """ The diagonal blocks of A, as the pass which built the system gathered
    them.

    Jacobi preconditioning, which is what Hendrickson and Konnert prescribe and
    what solve() falls back on, is close to useless on these systems, leaving
    the condition number within a factor of two of the unpreconditioned one;
    blocking a scatterer's parameters together instead takes it down by two to
    three orders of magnitude.

    max_block_size is not consulted here. The blocks were fixed before the pass
    that filled them, so this reports what was gathered rather than choosing
    afresh; the argument stays for the signature the callers share.
    """
    result, at = [], 0
    for indices in self.blocks:
      n = len(indices)
      block = self._blocks[at:at + n*n].reshape(n, n).copy()
      at += n*n
      if self.extra_design_matrix is not None:
        rows = self.extra_design_matrix[:, indices]
        block += rows.T.dot(
          self.extra_weights.as_numpy_array()[:, None]*rows)
      result.append((indices, block))
    return result

  def normal_matrix_diagonal(self):
    """ diag(A), for the Jacobi fallback in solve().

    Assembled from the blocks rather than from a pass of its own, which covers
    every parameter: parameter_blocks puts all of them in some block.
    """
    diagonal = flex.double(self.n_parameters, 0)
    at = 0
    for indices in self.blocks:
      n = len(indices)
      block = self._blocks[at:at + n*n].reshape(n, n)
      at += n*n
      for j, i in enumerate(indices):
        diagonal[int(i)] = block[j, j]
    if self.extra_design_matrix is not None:
      rows = self.extra_design_matrix
      extra = (self.extra_weights.as_numpy_array()[:, None]*rows*rows).sum(0)
      diagonal += flex.double(extra)
    return diagonal


class normal_matrix_problem(linearised_problem):
  """ The same system again, taken from the normal matrix rather than from J.

  The three ways of holding the Gauss-Newton system differ in what they store
  and in what they must then do with it, with m the number of conjugate
  gradient iterations:

      stored J      1 pass, no accumulation, CG over J   n_refl n_par
      normal matrix 1 pass, accumulation,    CG over A   n_par^2/2
      matrix free   m+1 passes                           n_par

  J is n_refl by n_par where A is n_par^2/2, so A is the smaller of the two
  until the parameters outnumber twice the reflections, which for ordinary data
  they never do. J also grows with the data rather than with the model, making
  it the first of the two to become impossible.

  Against Gauss-Newton, which stores the same matrix, this replaces the
  Cholesky decomposition with conjugate gradients: O(m n_par^2) rather than
  O(n_par^3/3). The decomposition is negligible for small structures, but it is
  the fastest growing term in the refinement and eventually the largest.

  The floating origin restraint is left in the build and project() deals with
  it as it does for the other modes, it contributing nothing to the projected
  operator, P (lambda d (x) d) P being zero when P d is. Leaving it out would
  cost two further passes over the reflections: one to learn the singular
  directions before the run and one to rebuild an invertible matrix for the
  covariance after it.

  It is not free, though. The restraint's weight lands on the diagonal of every
  parameter its directions touch -- for a polar axis, one coordinate of every
  atom in the structure -- so the blocks over-damp those directions and the
  conjugate gradients need several times as many iterations as they would with
  clean blocks. That is still the better trade, an iteration being orders of
  magnitude cheaper than a pass over the reflections, but the margin is a good
  deal narrower than the count of passes alone suggests. Capping
  max_cg_iterations bounds it for anyone who would rather not rely on that.
  """

  def __init__(self, ls, scale_factor, build=None):
    # build is how to make the normal equations. It matters that a caller can
    # supply it: normal_eqns_solving decorates the L.S. object to keep the
    # history of a run, and a build which goes straight to the undecorated one
    # is a cycle the history never hears about. This mode does the only builds
    # there are, so that would be the whole history.
    ls.overridden_scale_factor = scale_factor
    try:
      (build if build is not None else ls.build_up)()
    finally:
      ls.overridden_scale_factor = None

    self.n_parameters = ls.reparametrisation.n_independents
    self.packed_normal_matrix = \
      ls.normal_matrix_packed_u().deep_copy().as_numpy_array()
    self._right_hand_side = \
      ls.step_equations().right_hand_side().deep_copy()

    self.observables = ls.fc_sq.data()
    self.weights = ls.weights
    self.f_calc = ls.f_calc.data()
    self.scale_factor = ls.optimal_scale_factor()
    yo = ls.observations.fo_sq.data()
    self.sum_w_yo_sq = flex.sum(self.weights*yo*yo)
    self.weighted_observables_sq = flex.sum(
      self.weights*self.observables*self.observables)
    self.residuals = yo - self.scale_factor*self.observables

    # the restraints are already in A and b, build_up having added them, so
    # there are no extra rows to carry here
    self.extra_design_matrix = None
    self.extra_residuals = None
    self.extra_weights = None
    self.deflation_basis = None
    self._add_floating_origin(ls)

  def right_hand_side(self):
    return self._right_hand_side.deep_copy()

  def undamped_normal_matrix_times(self, v):
    """ A v, with A packed upper triangular as lstbx leaves it.

    Row major upper packing of a symmetric matrix is column major lower
    packing of the same one, which is what BLAS calls uplo='L'.
    """
    from scipy.linalg.blas import dspmv
    return flex.double(dspmv(n=self.n_parameters,
                             alpha=1.0,
                             ap=self.packed_normal_matrix,
                             x=v.as_numpy_array(),
                             lower=1))

  def _packed_index(self, i, j):
    if i > j: i, j = j, i
    return i*self.n_parameters - i*(i - 1)//2 + (j - i)

  def raw_blocks(self, ls, max_block_size=24):
    import numpy
    blocks = []
    for indices in self.parameter_blocks(ls, max_block_size):
      n = len(indices)
      block = numpy.empty((n, n))
      for a in range(n):
        for c in range(a, n):
          value = self.packed_normal_matrix[
            self._packed_index(int(indices[a]), int(indices[c]))]
          block[a, c] = block[c, a] = value
      blocks.append((indices, block))
    return blocks

  def normal_matrix_diagonal(self):
    diagonal = flex.double(self.n_parameters, 0)
    for i in range(self.n_parameters):
      diagonal[i] = self.packed_normal_matrix[self._packed_index(i, i)]
    return diagonal


class matrix_free_problem(linearised_problem):
  """ The same system, with the design matrix recomputed instead of stored.

  linearised_problem keeps the design matrix, which is n_refl by n_par and so
  outgrows the normal matrix it exists to avoid. Past that point the rows have
  to be recomputed inside each product, which is what ShelXL does.

  The arithmetic is identical -- the same builder, the same reflection loop,
  the same treatment of twinning, the solvent mask and the Fc corrections --
  and only the accumulator differs. Two of them live in
  smtbx/refinement/least_squares_matrix_free.h:

    - a summary pass, once per cycle, giving k*, grad k*, the right hand side
      and the diagonal blocks
    - a product pass, once per conjugate gradient iteration, giving A p

  Each holds O(n_par) per thread against the O(n_par^2) a normal matrix needs,
  which is what lets this run at sizes where forming A is out of the question.
  The cost is one pass over the reflections per iteration rather than one per
  cycle, which is a great many more, so this is a way to run where the others
  cannot rather than a way to run faster. See cgls_iterations.mode_for.
  """

  def __init__(self, ls, scale_factor, f_mask_data=None):
    from smtbx_refinement_least_squares_ext import (
      separable_scale_factor_summary)
    reparametrisation = ls.reparametrisation
    self.fc_correction = reparametrisation.fc_correction
    if self.fc_correction is None:
      self.fc_correction = xray.dummy_fc_correction()
    if f_mask_data is None:
      f_mask_data = ls.mask_data()
    self.ls = ls
    self.f_mask_data = f_mask_data
    self.weighting_scheme = ls.weighting_scheme
    self.weighting_scale_factor = scale_factor
    self.jacobian_transpose = \
      reparametrisation.jacobian_transpose_matching_grad_fc()
    self.n_parameters = self.jacobian_transpose.n_rows
    self.may_parallelise = getattr(ls, 'may_parallelise', False)

    # The blocks have to be settled before the pass which fills them, the
    # accumulator gathering each block's own outer product as it goes.
    self.blocks = self.parameter_blocks(ls)
    self.block_parameters = flex.int([int(i) for b in self.blocks for i in b])
    self.block_sizes = flex.int([len(b) for b in self.blocks])

    summary = separable_scale_factor_summary(
      self.n_parameters, self.block_parameters, self.block_sizes)
    built = self._build(summary)

    self.observables = flex.double(built.observables())
    self.weights = flex.double(built.weights())
    self.f_calc = built.f_calc()
    self.sum_w_yo_sq = summary.sum_w_yo_sq()
    self.weighted_observables_sq = summary.sum_w_yc_sq()
    self.scale_factor = summary.scale_factor()
    self.grad_scale_factor = summary.grad_scale_factor()
    self.residuals = (ls.observations.fo_sq.data()
                      - self.scale_factor*self.observables)
    self._right_hand_side = summary.right_hand_side()
    self._blocks = summary.blocks().as_numpy_array()

    self.extra_design_matrix = None
    self.extra_residuals = None
    self.extra_weights = None
    self.deflation_basis = None
    self._add_restraints(ls)
    self._add_floating_origin(ls)

  def _build(self, accumulator):
    """ One pass over the reflections, accumulating into the given object. """
    from smtbx_refinement_least_squares_ext import build_normal_equations
    return build_normal_equations(
      accumulator, self.ls.observations, self.f_mask_data,
      self.weighting_scheme, self.weighting_scale_factor,
      self.ls.one_h_linearisation, self.jacobian_transpose,
      self.fc_correction, False, self.may_parallelise, False, 300)

  # right_hand_side, raw_blocks and normal_matrix_diagonal are the inherited
  # ones: both modes get them from the same summary accumulator, the only
  # difference being that this one keeps no design matrix beside it.

  def undamped_normal_matrix_times(self, v):
    from smtbx_refinement_least_squares_ext import (
      separable_scale_factor_product)
    product = separable_scale_factor_product(
      self.n_parameters, self.scale_factor, self.grad_scale_factor, v,
      self.sum_w_yo_sq)
    self._build(product)
    result = product.result()
    if self.extra_design_matrix is not None:
      rows = self.extra_design_matrix
      u = rows.dot(v.as_numpy_array())*self.extra_weights.as_numpy_array()
      result += flex.double(rows.T.dot(u))
    return result


def report_state(problem, ls):
  """ Give the L.S. object what it needs to report a cycle.

  wR2, the goodness of fit and R1 are ordinarily computed from quantities
  build_up leaves behind. Nothing here builds the normal equations, so those
  quantities are handed over from the linearised problem instead, which has all
  of them. The alternative would be an extra pass over the reflections per
  cycle purely to report on it.
  """
  observations = ls.observations
  ls.f_calc = observations.fo_sq.array(data=problem.f_calc, sigmas=None)
  ls.fc_sq = observations.fo_sq.array(
    data=problem.observables, sigmas=None).set_observation_type_xray_intensity()
  ls.weights = problem.weights
  ls.overridden_scale_factor = problem.scale_factor

  weighted_residual = flex.sum(problem.weights*flex.pow2(problem.residuals))
  ls.objective_data_only = weighted_residual/(2*problem.sum_w_yo_sq)
  degrees_of_freedom = problem.residuals.size() - problem.n_parameters
  if degrees_of_freedom > 0:
    ls.chi_sq_data_only = weighted_residual/degrees_of_freedom
    restrained = weighted_residual
    if problem.extra_residuals is not None:
      restrained += flex.sum(
        problem.extra_weights*flex.pow2(problem.extra_residuals))
      degrees_of_freedom += problem.extra_residuals.size()
    ls.chi_sq_data_and_restraints = restrained/degrees_of_freedom


def solve(problem, max_iterations=None, tolerance=1e-4, preconditioner=None):
  """ Solve A x = b by preconditioned conjugate gradients.

  Returns the step and the number of iterations it took. Stopping on the
  residual relative to the right hand side is the usual criterion; the default
  is loose because the system is only a local quadratic model of the objective
  and solving it to great accuracy buys nothing.

  preconditioner is the list of blocks from block_preconditioner; passing None
  falls back on Jacobi, which is much the worse of the two.
  """
  b = problem.project(problem.right_hand_side())
  n = b.size()
  if max_iterations is None: max_iterations = n

  if preconditioner is None:
    diagonal = problem.normal_matrix_diagonal()
    diagonal.set_selected(diagonal <= 0, 1.)
    def apply_preconditioner(v): return v/diagonal
  else:
    def apply_preconditioner(v):
      a = v.as_numpy_array()
      out = flex.double(n, 0)
      result = out.as_numpy_array()
      for indices, inverse in preconditioner:
        result[indices] = inverse.dot(a[indices])
      return flex.double(result)

  x = flex.double(n, 0)
  r = b.deep_copy()
  z = problem.project(apply_preconditioner(r))
  p = z.deep_copy()
  r_dot_z = flex.sum(r*z)
  b_norm = b.norm()
  if b_norm == 0: return x, 0

  for iteration in range(1, max_iterations + 1):
    ap = problem.project(problem.normal_matrix_times(p))
    p_dot_ap = flex.sum(p*ap)
    if p_dot_ap <= 0:
      # A is positive definite in exact arithmetic; losing that means the
      # iteration has run out of numerical room and the step in hand is the
      # best available
      break
    alpha = r_dot_z/p_dot_ap
    x += alpha*p
    r -= alpha*ap
    if r.norm() <= tolerance*b_norm: return x, iteration
    z = problem.project(apply_preconditioner(r))
    r_dot_z_next = flex.sum(r*z)
    p *= r_dot_z_next/r_dot_z
    p += z
    r_dot_z = r_dot_z_next
  return x, iteration


class cgls_iterations(normal_eqns_solving.iterations):
  """ Refine by conjugate gradient least squares.

  One cycle is one linearisation of the problem -- one pass over the
  reflections computing structure factors and their derivatives -- followed by
  conjugate gradients on the resulting linear system. That is the shape ShelXL
  gives its CGLS instruction, and it is what makes a cycle here comparable with
  a Gauss-Newton cycle rather than with one evaluation of an objective.

  What is held between the linearisation and the solve is set by mode below,
  and decides both how much memory a run needs and where its time goes.

  Standard uncertainties cannot come out of the conjugate gradients themselves,
  there being no normal matrix to invert, so a run ends with normal equations
  standing for the covariance matrix to be taken from. ShelXL has the same
  limitation and the same remedy.
  """

  n_max_iterations = 20
  max_cg_iterations = None
  # How exactly to solve the linearised system each cycle. It is only a local
  # model of the objective, so solving it tighter than the model is worth buys
  # nothing, which is why ShelXL's CGLS solves it loosely and this default is
  # 1e-3. In stored mode a looser solve directly saves passes over the design
  # matrix, most of all on ill-conditioned systems. See notes/bench_tol.py.
  cg_tolerance = 1e-3
  preconditioner_block_size = 24
  # Below this the shifts are small enough to call it converged. Measured
  # against the parameters themselves rather than their s.u., which are not
  # available while the normal matrix is not being formed.
  convergence_as_relative_shift = 1e-4
  # The first argument of ShelXL's DAMP, on the scale lstbx uses for it, which
  # is a thousandth of the one the instruction is written on. The same default
  # normal_eqns_solving.naive_iterations_with_damping takes, so that the two
  # respond to the instruction identically.
  damping_value = 0.0007
  # How the linearised system is held: 'stored' keeps the design matrix,
  # 'normal_matrix' keeps A instead, 'matrix_free' keeps neither, and None
  # chooses on what will fit.
  #
  # Storing the design matrix is the quickest, being the only one of the three
  # which skips the rank-one accumulation, but it is also much the largest,
  # n_refl by n_par against the normal matrix's n_par^2/2, and it grows with the
  # data rather than with the model. So it is the first to become impossible,
  # and what replaces it keeps the same single pass over the reflections per
  # cycle while swapping the Cholesky decomposition for conjugate gradients.
  #
  # That swap is what makes this scale. The decomposition is a small part of a
  # cycle for a small structure, but it grows as O(n_par^3) where the reflection
  # loop grows as O(n_refl n_atoms), so it eventually dominates everything else.
  # Conjugate gradients make it O(m n_par^2).
  mode = None
  # Above this many MB a stored design matrix is judged not to fit. None asks
  # the machine instead, which is what this should mean: a fixed number is a
  # guess about hardware, and the mode it excludes is much the quickest one. The
  # rule has to mean "the design matrix will not fit", never "it is large".
  max_design_matrix_memory = None
  # What share of the machine a design matrix may take when the above is None. A
  # quarter leaves room for everything else resident, all of it far smaller than
  # the matrix at any size where this matters.
  design_matrix_memory_fraction = 0.25
  # And what to allow when the machine cannot be asked. Generous on purpose: a
  # too-large allocation announces itself, a too-shy budget silently picks a
  # slower mode.
  design_matrix_memory_if_unknown = 8192
  # Whether to end the run with the normal equations standing, the only way the
  # standard uncertainties can be had, since the conjugate gradients never form
  # the matrix they come from. That closing build is a full Gauss-Newton one, so
  # turn it off for a refinement whose result is another refinement, as SHELXL's
  # CGLS does.
  #
  # Named for what it does and *not* `standard_uncertainties`: that is a method
  # on the olex2 iteration classes this mixes into, and a keyword of that name
  # lands on the instance through adopt_optional_init_args and shadows the
  # method. A flag must never take the name of a method the same object carries.
  compute_standard_uncertainties = True
  # Hold the design matrix in single precision, accumulating the products
  # against it in double. The matrix is half the size and, the products being
  # bandwidth bound, so is the solve; the double accumulation keeps the error
  # negligible. Off by default because it moves the last digits of a result
  # someone may already have. It also halves what the matrix must fit within the
  # budget, so it can decide the mode.
  single_precision_design_matrix = False
  # And above this, neither does the normal matrix, at which point there is
  # nothing left but to recompute.
  max_normal_matrix_memory = 8192
  verbose = False
  log = None

  def design_matrix_budget(self):
    """ How many MB a stored design matrix may occupy.

    Asks the machine unless told a number. Where it cannot be asked, it returns
    the generous fallback rather than a small one: this budget gates the
    quickest mode, so a cautious guess here does not fail safely, it silently
    refuses it.
    """
    if self.max_design_matrix_memory is not None:
      return self.max_design_matrix_memory
    from libtbx.utils import guess_total_memory
    total = guess_total_memory()
    if total is None:
      return self.design_matrix_memory_if_unknown
    return total/1048576.*self.design_matrix_memory_fraction

  # Why mode_for chose what it chose, for a caller which reports it. A mode
  # picked automatically and reported without its reason is how a structure
  # went on taking the slow path unnoticed.
  mode_reason = None

  def mode_for(self, ls):
    """ Which of the three ways to hold the system, given its size. """
    single = self.single_precision_design_matrix
    precision = "single precision" if single else "double precision"
    if self.mode is not None:
      # The precision is worth saying even where the mode was not chosen here.
      # It sets the size of the matrix and moves the last digits of the result,
      # so a log which reports the mode and not the precision leaves out half
      # of what distinguishes two runs -- which is the mistake the budget made
      # before it began reporting its reason at all.
      self.mode_reason = ("asked for, %s" % precision if self.mode == 'stored'
                          else "asked for")
      return self.mode
    n_refl = ls.observations.fo_sq.size()
    n_par = ls.reparametrisation.n_independents
    element = linearised_problem.bytes_per_element[bool(single)]
    wanted = n_refl*n_par*element/1048576.
    budget = self.design_matrix_budget()
    self.mode_reason = ("design matrix %.0f MB (%s) against a budget of %.0f MB"
                        % (wanted, precision, budget))
    if wanted <= budget:
      return 'stored'
    if n_par*(n_par + 1)/2*8/1048576. <= self.max_normal_matrix_memory:
      return 'normal_matrix'
    return 'matrix_free'

  def problem_for(self, ls, scale_factor):
    mode = self.mode_for(ls)
    # kept so that a caller can report which of the three actually ran, auto
    # being free to choose differently from one structure to the next
    self.chosen_mode = mode
    if mode == 'stored':
      return linearised_problem(
        ls, scale_factor,
        single_precision=self.single_precision_design_matrix)
    if mode == 'normal_matrix':
      return normal_matrix_problem(ls, scale_factor,
                                   build=self.non_linear_ls.build_up)
    if mode == 'matrix_free':
      return matrix_free_problem(ls, scale_factor)
    raise RuntimeError("unrecognised mode %r" % (mode,))

  def do(self):
    ls = self.non_linear_ls.actual
    self.n_iterations = 0
    self.n_cg_iterations = 0
    # Two builds bracket the run: one to begin, which settles the normalisation
    # of the restraints and makes the floating origin directions known, and one
    # to end, which leaves the normal equations standing for the covariance
    # matrix that conjugate gradients cannot produce.
    #
    # The normal matrix mode needs neither, building the equations itself once
    # per cycle and in full: the first cycle's build settles what the opening
    # one would have, and the last cycle's is still standing when the run ends.
    # That saves two passes over the reflections out of six, which is worth
    # having when a pass is the dominant cost of a cycle.
    bracketed = self.mode_for(ls) != 'normal_matrix'
    if bracketed:
      self.open_the_run(ls)
      scale_factor = ls.scale_factor()
    else:
      # the first build works one out for itself, having none to be told
      scale_factor = None

    while self.n_iterations < self.n_max_iterations:
      problem = self.problem_for(ls, scale_factor)
      problem.set_damping(self.damping_value)
      preconditioner = problem.block_preconditioner(
        ls, self.preconditioner_block_size)
      step, cg_iterations = solve(problem,
                                  max_iterations=self.max_cg_iterations,
                                  tolerance=self.cg_tolerance,
                                  preconditioner=preconditioner)
      self.n_cg_iterations += cg_iterations
      report_state(problem, ls)
      self.latest_problem = problem
      self.latest_preconditioner = preconditioner
      self.latest_step = step
      # the second argument of ShelXL's DAMP, and it has to happen before the
      # structure moves rather than after: a step is limited or it is not, and
      # reporting on one already taken is too late to limit it
      step = self.scale_shifts(step, problem, ls)
      self.latest_step = step

      self.begin_cycle()
      ls.apply_shifts(step)
      self.n_iterations += 1
      self.end_cycle()
      # what the weights of the next linearisation are computed with, exactly
      # as build_up would use the scale factor of the cycle before
      scale_factor = problem.scale_factor
      # after end_cycle, so that a caller reporting on the cycle has had its
      # say and can base the decision on what it worked out
      if self.had_small_enough_shifts(step, ls): break

    # A run has to close on a pass over the final model, because report_state
    # describes the model each cycle *entered* with -- without a closing pass
    # the objective and structure factors would be one shift out of date.
    # Standard uncertainties need the normal equations, so with them the close
    # is a full Gauss-Newton build; without them a gradient-free pass suffices.
    # The scale-factor override stays in place either way, so the final wR2 and
    # goodness of fit are weighted with the scale factor the last cycle used.
    if bracketed:
      self.non_linear_ls.build_up(
        objective_only=not self.compute_standard_uncertainties)
    ls.overridden_scale_factor = None

  def open_the_run(self, ls):
    """ Settle what the first cycle needs and cannot work out for itself.

    Two things: the factor which puts the restraints on the scale of the
    average residual, which wants the objective over the data alone, and the
    directions along which the structure may float, which project() strips out
    of every step.

    Neither wants a normal matrix. The objective comes from a gradient-free
    pass, and the floating origin directions are a function of the space group,
    the scatterers and the reparametrisation, not of the data -- the normal
    equations origin_fixing.add_to writes into only receive the restraint
    equation this mode discards in favour of projecting, so it is handed a
    scratch object and what it writes there is dropped. That leaves one
    gradient-free pass where a full Gauss-Newton build would otherwise run.
    """
    from scitbx.lstbx import normal_eqns
    self.non_linear_ls.build_up(objective_only=True)
    reparametrisation = ls.reparametrisation
    # add_to returns at once where the origin is fixed by the symmetry, which
    # is the common case; the scratch equations cost n_par^2/2 until then
    ls.origin_fixing_restraint.add_to(
      normal_eqns.linear_ls(reparametrisation.n_independents),
      reparametrisation.jacobian_transpose_matching_grad_fc(),
      reparametrisation.asu_scatterer_parameters)

  def scale_shifts(self, step, problem, ls):
    """ The second argument of ShelXL's DAMP: hold the largest shift/su to it.

    A no-op here, there being no standard uncertainties to measure a shift
    against while the normal matrix is not being formed. A caller which has an
    estimate of them, from the preconditioner's blocks or otherwise, overrides
    this and applies the limit as the Gauss-Newton path does.
    """
    return step

  def had_small_enough_shifts(self, step, ls):
    norm = ls.reparametrisation.norm_of_independent_parameter_vector
    return step.norm() <= self.convergence_as_relative_shift*max(norm, 1e-30)

  def begin_cycle(self):
    """ Hook for a caller which reports on each cycle; a no-op here. """

  def end_cycle(self):
    """ Hook called once the shifts of a cycle have been applied. """
    if not self.verbose: return
    out = self.log if self.log is not None else sys.stdout
    ls = self.non_linear_ls.actual
    print("  %4i  %12.8f  %8.4f  %8.4f  %6d"
          % (self.n_iterations, ls.objective_data_only, ls.wR2(),
             ls.r1_factor()[0], self.n_cg_iterations), file=out)

  def __str__(self):
    mode = getattr(self, 'chosen_mode', None)
    if mode is None:
      return "conjugate gradient least squares"
    if self.mode_reason is None:
      return "conjugate gradient least squares, %s" % mode
    return ("conjugate gradient least squares, %s (%s)"
            % (mode, self.mode_reason))
