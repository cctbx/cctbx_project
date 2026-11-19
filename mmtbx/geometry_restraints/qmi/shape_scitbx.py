"""Coordination-geometry classification by Continuous Shape Measures (plan 5C).

For a metal's first coordination shell, assign the closest ideal polyhedron and
report how good the fit is, feeding the electronic-state guess (plan 6). The
measure is the Continuous Shape Measure S(Q,P) [11] (the polyhedral form of the
continuous symmetry measure [31]; minimal-distortion pathways [32]; shape maps
[33]): the residual after the best joint optimisation over rotation, isotropic
scale, and vertex permutation mapping the observed shell onto an ideal
polyhedron. S in [0, 100], 0 = perfect, larger = more distorted.

Two design choices specific to qmi:

- **Angular-only.** The observed donor directions are unit-normalised before the
  Procrustes, so S depends only on the *directions* to the donors, not the bond
  lengths. A Jahn-Teller-elongated or oxidation-shortened bond therefore cannot
  shift the geometry class, which is what keeps geometry (5C) separate from the
  state guess (6). The runtime path is metal-centred and origin-normalised; a
  ``centroid_normalized`` flag reproduces the canonical Pinsky-Avnir convention
  used to validate against published S values.

- **Carboxylate regularisation.** A chelating Asp/Glu (or a main-chain
  C-terminal carboxylate) donates two O to the same metal at a ~55-65 degree
  bite angle, which matches no regular-polyhedron edge and would force a
  spurious high-CN misassignment. A symmetric chelate (two near-equal M-O) is
  collapsed to one effective vertex along the O-O bisector for the shape fit
  only; an anisobidentate carboxylate (unequal M-O) keeps only its near O, the
  far O being a long contact. The unweighted bisector is deliberate: distance
  weighting would reintroduce the bond-length dependence the angular-only
  design excludes. ``cn`` (effective) is for geometry only; the electronic-state
  guess (plan 6) uses the coordinative count ``cn + n_bidentate_carboxylate`` (a
  symmetric chelate restored to its two O, an anisobidentate far O left out as a
  long contact), while ``raw_cn`` (every shell donor) is the emitted M-L bond
  count. The QM charge/cluster still sees both O.

Coordinate math uses only scitbx flex arrays and scitbx.linalg (no numpy). This
is the pure-scitbx parallel of shape.py (the numpy-batched version used in the
pipeline); the two produce identical results and are kept side by side for
comparison. The permutation search is a Python loop with one 3x3 scitbx SVD
(singular values only) per permutation, over the symmetry-reduced set by
default.
"""

from math import sqrt, cos, sin, pi
from itertools import permutations
from functools import lru_cache

from cctbx.array_family import flex
from scitbx import matrix
from scitbx.linalg import svd
from libtbx import group_args

from mmtbx.geometry_restraints.qmi import util


# ---------------------------------------------------------------------------
# Ideal polyhedra: unit-radius vertices centred at the origin. These are
# mathematics (analytic vertex sets), so they live in code, not a data file.
# ---------------------------------------------------------------------------

def _pentagon_xy(z):
  """Five vertices of a regular pentagon in the plane z, each at unit length."""
  r = sqrt(1.0 - z * z)
  return [(r * cos(2.0 * pi * k / 5.0), r * sin(2.0 * pi * k / 5.0), z)
          for k in range(5)]


def _ideal_linear():
  return flex.vec3_double([(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)])


def _ideal_trigonal_planar():
  return flex.vec3_double([
    (1.0, 0.0, 0.0),
    (-0.5, sqrt(3.0) / 2.0, 0.0),
    (-0.5, -sqrt(3.0) / 2.0, 0.0)])


def _ideal_tetrahedral():
  s = 1.0 / sqrt(3.0)
  return flex.vec3_double([(s, s, s), (s, -s, -s), (-s, s, -s), (-s, -s, s)])


def _ideal_square_planar():
  return flex.vec3_double([
    (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0)])


def _ideal_trigonal_bipyramidal():
  return flex.vec3_double([
    (1.0, 0.0, 0.0),
    (-0.5, sqrt(3.0) / 2.0, 0.0),
    (-0.5, -sqrt(3.0) / 2.0, 0.0),
    (0.0, 0.0, 1.0), (0.0, 0.0, -1.0)])


def _ideal_square_pyramidal():
  # Four basal vertices (z=0) plus an apex on +z: the vacant-octahedron
  # idealisation (90 degree apex-metal-basal angle), a common square-pyramid
  # reference. SHAPE's SPY-5 opens the apex further; not used here, as no
  # verified apex angle was available.
  return flex.vec3_double([
    (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0), (0.0, 1.0, 0.0), (0.0, -1.0, 0.0),
    (0.0, 0.0, 1.0)])


def _ideal_octahedral():
  return flex.vec3_double([
    (1.0, 0.0, 0.0), (-1.0, 0.0, 0.0),
    (0.0, 1.0, 0.0), (0.0, -1.0, 0.0),
    (0.0, 0.0, 1.0), (0.0, 0.0, -1.0)])


def _ideal_trigonal_prism():
  # Equilateral (all-edges-equal) D3h prism, the canonical SHAPE TPR-6:
  # r = 2/sqrt(7), h = sqrt(3/7), r^2 + h^2 = 1.
  r = 2.0 / sqrt(7.0)
  h = sqrt(3.0 / 7.0)
  top = [(r * cos(2.0 * pi * k / 3.0), r * sin(2.0 * pi * k / 3.0), h)
         for k in range(3)]
  bot = [(x, y, -h) for (x, y, _) in top]
  return flex.vec3_double(top + bot)


def _ideal_pentagonal_bipyramidal():
  return flex.vec3_double(
    _pentagon_xy(0.0) + [(0.0, 0.0, 1.0), (0.0, 0.0, -1.0)])


def _ideal_square_antiprism():
  # Uniform (all-edges-equal) D4d antiprism, the canonical SHAPE SAPR-8:
  # R = sqrt(4/(4+sqrt(2))), h = sqrt(1-R^2); vertices unit length.
  r = sqrt(4.0 / (4.0 + sqrt(2.0)))
  h = sqrt(1.0 - r * r)
  top = [(r * cos(pi / 2.0 * k), r * sin(pi / 2.0 * k), h) for k in range(4)]
  bot = [(r * cos(pi / 4.0 + pi / 2.0 * k),
          r * sin(pi / 4.0 + pi / 2.0 * k), -h) for k in range(4)]
  return flex.vec3_double(top + bot)


def _ideal_cubic():
  s = 1.0 / sqrt(3.0)
  return flex.vec3_double([
    (s, s, s), (s, s, -s), (s, -s, s), (s, -s, -s),
    (-s, s, s), (-s, s, -s), (-s, -s, s), (-s, -s, -s)])


IDEAL_POLYHEDRA = {
  util.LINEAR: _ideal_linear(),
  util.TRIGONAL: _ideal_trigonal_planar(),
  util.TETRAHEDRAL: _ideal_tetrahedral(),
  util.SQUARE_PLANAR: _ideal_square_planar(),
  util.TRIGONAL_BIPYRAMIDAL: _ideal_trigonal_bipyramidal(),
  util.SQUARE_PYRAMIDAL: _ideal_square_pyramidal(),
  util.OCTAHEDRAL: _ideal_octahedral(),
  util.TRIGONAL_PRISM: _ideal_trigonal_prism(),
  util.PENTAGONAL_BIPYRAMIDAL: _ideal_pentagonal_bipyramidal(),
  util.SQUARE_ANTIPRISM: _ideal_square_antiprism(),
  util.CUBIC: _ideal_cubic(),
}


# Candidate shapes worth scoring at each coordination number.
SHAPES_BY_CN = {
  2: [util.LINEAR],
  3: [util.TRIGONAL],
  4: [util.TETRAHEDRAL, util.SQUARE_PLANAR],
  5: [util.TRIGONAL_BIPYRAMIDAL, util.SQUARE_PYRAMIDAL],
  6: [util.OCTAHEDRAL, util.TRIGONAL_PRISM],
  7: [util.PENTAGONAL_BIPYRAMIDAL],
  8: [util.SQUARE_ANTIPRISM, util.CUBIC],
}

# Brute-force permutation search is over N! labelings; 8! = 40320 is trivial
# vectorised, 9! is large. CN above this is flagged, not classified.
_MAX_BRUTE_N = 8

# Flag a shape assignment as ambiguous when the two best fits are within this
# many CShM units (heuristic; distorted intermediates sit between polyhedra).
AMBIGUITY_GAP = 2.0


# ---------------------------------------------------------------------------
# Symmetry-reduced permutations (the default; pass symmetry_reduced=False for
# the full N! search). The CShM over all N! vertex labelings is redundant: if g
# is a proper-rotation symmetry of the ideal polyhedron (P[g] = R P, det R = +1)
# then labelings pi and g.pi score identically (the rotation R is absorbed by
# the Procrustes), so only one labeling per orbit of the automorphism group need
# be scored. Reductions e.g. octahedron 720 -> 30, cube 40320 -> 1680.
# Precomputed once per shape from the analytic vertices (pure Python). The
# centroid-normalised validation path always uses the full N! set.
# ---------------------------------------------------------------------------

_SYM_TOL = 1e-6


def _dot3(a, b):
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]


def _det3(a, b, c):
  return (a[0] * (b[1] * c[2] - b[2] * c[1])
          - a[1] * (b[0] * c[2] - b[2] * c[0])
          + a[2] * (b[0] * c[1] - b[1] * c[0]))


def _gram_preserving_perms(V):
  """All permutations g preserving every pairwise dot product of the vertex
  set V (its full orthogonal automorphism group), by constraint backtracking."""
  n = len(V)
  gram = [[_dot3(V[i], V[j]) for j in range(n)] for i in range(n)]
  out = []
  assign = [-1] * n
  used = [False] * n

  def bt(k):
    if k == n:
      out.append(tuple(assign))
      return
    for j in range(n):
      if used[j]:
        continue
      if all(abs(gram[i][k] - _dot3(V[assign[i]], V[j])) <= _SYM_TOL
             for i in range(k)):
        assign[k] = j
        used[j] = True
        bt(k + 1)
        used[j] = False
        assign[k] = -1

  bt(0)
  return out


def _proper_symmetry_perms(V):
  """Gram-preserving permutations realisable by a proper rotation. For a
  planar/linear vertex set all are proper (compose with the in-plane
  reflection, which fixes the vertices); for a 3D set, keep those preserving
  the sign of a non-coplanar triple's determinant."""
  perms = _gram_preserving_perms(V)
  n = len(V)
  tri = None
  for i in range(n):
    for j in range(i + 1, n):
      for k in range(j + 1, n):
        if abs(_det3(V[i], V[j], V[k])) > _SYM_TOL:
          tri = (i, j, k)
          break
      if tri:
        break
    if tri:
      break
  if tri is None:
    return perms
  i0, i1, i2 = tri
  base_pos = _det3(V[i0], V[i1], V[i2]) > 0
  return [g for g in perms
          if (_det3(V[g[i0]], V[g[i1]], V[g[i2]]) > 0) == base_pos]


def _orbit_representatives(n, group):
  """One permutation per orbit of S_n under left composition g.perm."""
  seen = set()
  reps = []
  for perm in permutations(range(n)):
    if perm in seen:
      continue
    reps.append(perm)
    for g in group:
      seen.add(tuple(g[perm[k]] for k in range(n)))
  return reps


@lru_cache(maxsize=None)
def _reduced_perms(shape_name):
  """Symmetry-reduced vertex labelings for a shape (see the section note)."""
  V = [tuple(v) for v in IDEAL_POLYHEDRA[shape_name]]
  return _orbit_representatives(len(V), _proper_symmetry_perms(V))


@lru_cache(maxsize=None)
def _all_perms(n):
  """All N! permutations of range(n) as a tuple of index tuples."""
  return tuple(permutations(range(n)))


# ---------------------------------------------------------------------------
# Continuous Shape Measure: scaled orthogonal Procrustes per vertex
# permutation, via scitbx/flex (no numpy), over the selected permutation set
# (full N! by default, or the symmetry-reduced set). One 3x3 scitbx SVD
# (singular values only) per permutation, in a Python loop.
# ---------------------------------------------------------------------------

def _covariance(P, Q):
  """3x3 covariance M[i][j] = sum_k P[k][i] * Q[k][j] as a scitbx matrix.sqr.

  Uses flex.vec3_double.parts() so the N-length sums run in C++.
  """
  pc = P.parts()
  qc = Q.parts()
  return matrix.sqr(tuple(
    flex.sum(pc[i] * qc[j]) for i in range(3) for j in range(3)))


def _csm_for_permutation(Q, P_pi, q_sq, p_sq):
  """CShM (0-100) for one fixed permutation of the ideal vertices.

  Scaled orthogonal Procrustes: singular values of M = P_pi^T Q give
  S = 100 (1 - (sum sigma)^2 / (|Q|^2 |P|^2)); the smallest sigma is
  sign-flipped when det(M) < 0 (a reflection) so the rotation is proper.
  Only singular values are needed, so u/v are not accumulated, and the
  reflection sign comes from det(M) directly (M = U S Vt, S >= 0).
  """
  M = _covariance(P_pi, Q)
  mf = flex.double(M.elems)
  mf.reshape(flex.grid(3, 3))
  sigma = list(svd.real(mf, accumulate_u=False, accumulate_v=False).sigma)
  smallest = -sigma[2] if M.determinant() < 0.0 else sigma[2]
  trace_max = sigma[0] + sigma[1] + smallest
  s = 100.0 * (1.0 - (trace_max * trace_max) / (q_sq * p_sq))
  return 0.0 if s < 0.0 else s


def _best_csm_over_permutations(Q, P, perms):
  """Minimum CShM of observed ``Q`` against ideal ``P`` over the given vertex
  permutations (an iterable of index tuples), both flex.vec3_double.
  Returns (best_csm, best_permutation_tuple)."""
  q_sq = flex.sum(Q.dot(Q))
  p_sq = flex.sum(P.dot(P))
  best_csm = None
  best_perm = None
  for perm in perms:
    P_pi = P.select(flex.size_t(perm))
    s = _csm_for_permutation(Q, P_pi, q_sq, p_sq)
    if best_csm is None or s < best_csm:
      best_csm = s
      best_perm = perm
  return best_csm, best_perm


def _centered(v):
  """flex.vec3_double centred on its own centroid."""
  return v - flex.vec3_double(v.size(), v.mean())


def csm_against_shape(observed, shape_name, centroid_normalized=False,
                      symmetry_reduced=True):
  """CShM of ``observed`` (flex.vec3_double) against a named ideal polyhedron.

  ``observed`` is metal-centred; qmi passes unit-normalised directions (see the
  module docstring). ``centroid_normalized=False`` (default) normalises around
  the origin (the metal); ``True`` is the canonical Pinsky-Avnir convention
  (centre both point sets on their centroids) used to compare against published
  S values.

  ``symmetry_reduced`` (default True) scores only one vertex labeling per
  automorphism-group orbit instead of all N! (identical result, fewer
  permutations); pass False for the full search. It is ignored on the
  centroid-normalised path, which always uses the full set.

  Returns group_args(shape, csm, gRMSD_deg, permutation). Raises ValueError on
  a vertex-count mismatch or N above the brute-force cap.
  """
  ideal = IDEAL_POLYHEDRA[shape_name]
  n = observed.size()
  if n != ideal.size():
    raise ValueError(f"shape {shape_name!r} has {ideal.size()} vertices but "
                     f"observed has {n}")
  if n > _MAX_BRUTE_N:
    raise ValueError(f"CN {n} exceeds the brute-force CShM cap "
                     f"({_MAX_BRUTE_N})")

  if centroid_normalized:
    Q = _centered(observed)
    P = _centered(ideal)
  else:
    Q = observed
    P = ideal
  if flex.sum(Q.dot(Q)) == 0.0:
    raise ValueError("observed coordinates are all zero after centering")

  perms = (_reduced_perms(shape_name)
           if symmetry_reduced and not centroid_normalized else _all_perms(n))
  csm, perm = _best_csm_over_permutations(Q, P, perms)
  # gRMSD on the original (uncentered) observed vs the best-permuted ideal;
  # angles are translation/scale invariant so the convention does not matter.
  P_best = ideal.select(flex.size_t(perm))
  grmsd = _angular_grmsd(observed, P_best)
  return group_args(shape=shape_name, csm=csm, gRMSD_deg=grmsd,
                    permutation=perm)


def find_best_shape(observed, candidate_shapes=None, centroid_normalized=False,
                    symmetry_reduced=True):
  """Score ``observed`` against every candidate shape for its N and return the
  best fit as group_args(shape, csm, gRMSD_deg, permutation, all_fits).

  ``candidate_shapes`` defaults to SHAPES_BY_CN[N]. ``symmetry_reduced`` is
  forwarded to :func:`csm_against_shape`. Returns None if N exceeds the
  brute-force cap or no candidate shape exists for N.
  """
  n = observed.size()
  if n > _MAX_BRUTE_N:
    return None
  if candidate_shapes is None:
    candidate_shapes = SHAPES_BY_CN.get(n)
    if not candidate_shapes:
      return None

  all_fits = {}
  best = None
  for shape in candidate_shapes:
    fit = csm_against_shape(observed, shape,
                            centroid_normalized=centroid_normalized,
                            symmetry_reduced=symmetry_reduced)
    all_fits[shape] = fit
    if best is None or fit.csm < best.csm:
      best = fit
  return group_args(shape=best.shape, csm=best.csm, gRMSD_deg=best.gRMSD_deg,
                    permutation=best.permutation, all_fits=all_fits)


# ---------------------------------------------------------------------------
# gRMSD: RMSD over sorted pairwise L-M-L angles (Zheng et al. 2013, ref [10])
# ---------------------------------------------------------------------------

def _pairwise_angles_deg(v):
  out = []
  n = v.size()
  for i in range(n):
    vi = matrix.col(v[i])
    for j in range(i + 1, n):
      out.append(vi.angle(matrix.col(v[j]), deg=True))
  out.sort()
  return out


def _angular_grmsd(observed, ideal_permuted):
  ao = _pairwise_angles_deg(observed)
  ai = _pairwise_angles_deg(ideal_permuted)
  if not ao:
    return 0.0
  ss = sum((ao[i] - ai[i]) ** 2 for i in range(len(ao)))
  return sqrt(ss / len(ao))


# ---------------------------------------------------------------------------
# Shell preparation: unit directions + carboxylate regularisation
# ---------------------------------------------------------------------------

# Carboxylate O pairs by residue; both O coordinating the same metal at a small
# bite is a chelating bidentate (see module docstring). Acetate (ACT: O/OXT) is
# already covered by the C-terminal O/OXT rule. Bidentate oxoanions with three+
# O (carbonate, phosphate/sulfate) need any-two-of-N logic and are future scope.
_CARBOXYLATE_O_PAIRS = {
  "ASP": ("OD1", "OD2"),
  "GLU": ("OE1", "OE2"),
  "FMT": ("O1", "O2"),    # formate ligand
}
# Main-chain C-terminal carboxylate (also acetate), for any residue name.
_CTERM_O_PAIR = ("O", "OXT")

# Above this O-M-O angle a "carboxylate" O pair is not a tight chelate; keep the
# O separate and flag rather than collapsing (defensive; syn chelation is ~60).
_CHELATE_MAX_ANGLE_DEG = 90.0

# A chelate has two near-equal M-O bonds. If the two M-O differ by more than
# this, the carboxylate is anisobidentate/monodentate: only the near O is a
# coordination vertex and the far O is a long contact. Boundary from CSD
# carboxylate binding-mode surveys (Rardin, Tolman, Lippard, New J. Chem. 1991,
# 15, 417).
_CHELATE_MAX_ASYMMETRY_A = 0.5


def _unit(vec):
  m = matrix.col(vec)
  n = m.length()
  return None if n == 0.0 else (m / n).elems


def _effective_donors(atoms, metal_xyz, donor_iseqs):
  """Build the effective donor unit-direction set for shape scoring.

  Groups donor O by carboxylate identity; a carboxylate with both O in the
  shell is collapsed to one bisector vertex. Returns
  group_args(units, cn, raw_cn, n_bidentate_carboxylate, flags).
  """
  m = matrix.col(metal_xyz)
  # Per-residue map: {carboxylate-O-name: iseq} for residues that could chelate.
  by_res = {}
  for iseq in donor_iseqs:
    a = atoms[iseq]
    resname = a.parent().resname.strip()
    name = a.name.strip()
    pair = _CARBOXYLATE_O_PAIRS.get(resname)
    if pair is not None and name in pair:
      by_res.setdefault(util.residue_key(a), {})[name] = iseq
    elif name in _CTERM_O_PAIR:
      by_res.setdefault(util.residue_key(a), {})[name] = iseq

  chelating = {}  # iseq -> bisector partner iseq (symmetric chelate)
  drop = set()    # far O of an anisobidentate carboxylate (a long contact)
  flags = []
  n_bidentate = 0
  for key, names in by_res.items():
    resname = key[-1]
    pair = _CARBOXYLATE_O_PAIRS.get(resname, _CTERM_O_PAIR)
    if pair[0] not in names or pair[1] not in names:
      continue
    i0, i1 = names[pair[0]], names[pair[1]]
    v0 = matrix.col(atoms[i0].xyz) - m
    v1 = matrix.col(atoms[i1].xyz) - m
    d0, d1 = v0.length(), v1.length()
    if d0 == 0.0 or d1 == 0.0:
      continue
    if v0.angle(v1, deg=True) > _CHELATE_MAX_ANGLE_DEG:
      flags.append("carboxylate_wide_bite")   # not a tight chelate; keep both
      continue
    if abs(d0 - d1) > _CHELATE_MAX_ASYMMETRY_A:
      drop.add(i0 if d0 > d1 else i1)          # far O is a long contact
      flags.append("carboxylate_asymmetric")
      continue
    chelating[i0] = i1
    chelating[i1] = i0
    n_bidentate += 1

  units = []
  done = set()
  for iseq in donor_iseqs:
    if iseq in done:
      continue
    if iseq in drop:
      done.add(iseq)
      continue
    if iseq in chelating:
      partner = chelating[iseq]
      u0 = _unit(matrix.col(atoms[iseq].xyz) - m)
      u1 = _unit(matrix.col(atoms[partner].xyz) - m)
      if u0 is not None and u1 is not None:
        bis = _unit(matrix.col(u0) + matrix.col(u1))
        if bis is not None:
          units.append(bis)
      done.add(iseq)
      done.add(partner)
    else:
      u = _unit(matrix.col(atoms[iseq].xyz) - m)
      if u is not None:
        units.append(u)
      done.add(iseq)

  if n_bidentate:
    flags.append("bidentate_carboxylate")
  return group_args(
    units=flex.vec3_double(units),
    cn=len(units), raw_cn=len(donor_iseqs),
    n_bidentate_carboxylate=n_bidentate, flags=flags)


def classify_geometry(atoms, metal_iseq, donor_iseqs, centroid_normalized=False,
                      symmetry_reduced=True):
  """Assign a coordination polyhedron to a metal's first shell (plan 5C).

  ``atoms`` is the sub-model's hierarchy atom array (observed input geometry);
  ``donor_iseqs`` are the first-shell donors from
  ``restraints.detect_first_shell``. Directions are unit-normalised
  (angular-only) and chelating carboxylates are regularised to one vertex
  before scoring. ``symmetry_reduced`` is forwarded to :func:`find_best_shape`.

  Returns group_args with: ``shape`` (name or None), ``csm``, ``gRMSD_deg``,
  ``permutation``, ``ambiguous`` (bool), ``all_fits`` (dict), ``cn`` (effective,
  after carboxylate regularisation), ``raw_cn``, ``n_bidentate_carboxylate``,
  and ``flags``.
  """
  shell = _effective_donors(atoms, atoms[metal_iseq].xyz, donor_iseqs)
  base = group_args(
    shape=None, csm=None, gRMSD_deg=None, permutation=None, ambiguous=False,
    all_fits={}, cn=shell.cn, raw_cn=shell.raw_cn,
    n_bidentate_carboxylate=shell.n_bidentate_carboxylate, flags=shell.flags)

  if shell.cn < 2:
    base.flags = base.flags + ["cn_too_low"]
    return base
  if shell.cn > _MAX_BRUTE_N or shell.cn not in SHAPES_BY_CN:
    base.flags = base.flags + ["unsupported_cn"]
    return base

  best = find_best_shape(shell.units, centroid_normalized=centroid_normalized,
                         symmetry_reduced=symmetry_reduced)
  if best is None:
    base.flags = base.flags + ["unsupported_cn"]
    return base

  ordered = sorted(best.all_fits.values(), key=lambda f: f.csm)
  base.shape = best.shape
  base.csm = best.csm
  base.gRMSD_deg = best.gRMSD_deg
  base.permutation = best.permutation
  base.all_fits = best.all_fits
  base.ambiguous = (len(ordered) > 1
                    and (ordered[1].csm - ordered[0].csm) < AMBIGUITY_GAP)
  return base
