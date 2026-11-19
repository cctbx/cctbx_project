"""Modified-Seminario force constants -> restraint sigmas (plan Sec 5F).

Reads the Cartesian Hessian xTB writes (``--hess``, a ``<namespace>.hessian``
file, Hartree/Bohr^2, full 3N x 3N) at the optimised geometry and derives a
per-restraint force constant by the modified-Seminario method (Allen, Payne &
Cole 2018 [65], improving Seminario 1996):

- **Bond M-L**: project the (symmetrised) inter-atomic Hessian block onto the
  bond unit vector, summed over the block's eigenpairs.
- **Angle L-M-L**: project each of the two metal-donor blocks onto the
  in-plane perpendicular, then combine with the bond lengths.

The force constant is mapped to a restraint sigma by harmonic/thermal
equivalence, ``sigma = sqrt(k_B T / k)``, the positional spread of a harmonic
well of curvature ``k`` at temperature ``T`` (so a soft axial bond gives a
looser sigma, as the plan wants).  Sigmas are clamped to defensible bounds;
pathological (<= 0) curvature returns the loose cap, never a tight false
restraint.

Only the metal-donor and donor-metal-donor blocks are used (the "partial"
inner-sphere Hessian of plan Sec 5F): the metal and its donors are free during
the constrained optimisation, so their gradient is ~0 at convergence and those
blocks are at a true local minimum, independent of the frozen scaffold.
"""

import math

from libtbx import group_args
from libtbx.utils import null_out
from scitbx.array_family import flex
from scitbx import matrix
from scitbx.linalg import eigensystem


# Unit conversions.
HARTREE_TO_KCAL = 627.5094740631
BOHR_TO_ANG = 0.52917721067
# Hessian force constant Hartree/Bohr^2 -> kcal/mol/Angstrom^2.
HARTREE_BOHR2_TO_KCAL_A2 = HARTREE_TO_KCAL / (BOHR_TO_ANG ** 2)
# Boltzmann constant, kcal/(mol K).
KB_KCAL = 0.0019872041
# sigma = sqrt(k_B T / k) is the thermal spread at the crystal's data-collection
# temperature; almost all X-ray (cryostream) and cryo-EM (vitrified) structures
# are ~100 K, so that is the default (the qm.temperature PHIL overrides it for a
# room-temperature dataset).
DEFAULT_TEMPERATURE_K = 100.0

# Sigma clamps (plan: sigmas are floors to widen; never emit a false-tight
# restraint from numerical noise).
DEFAULT_BOND_SIGMA_MIN = 0.005    # A
DEFAULT_BOND_SIGMA_MAX = 0.20     # A
DEFAULT_ANGLE_SIGMA_MIN = 0.5     # deg
DEFAULT_ANGLE_SIGMA_MAX = 15.0    # deg

# Mass-weighted Cartesian Hessian eigenvalue (omega^2, in Hartree/(Bohr^2 amu))
# -> harmonic wavenumber (cm^-1): nu = sqrt(omega^2) / (2 pi c). Derived from
# CODATA constants so the factor is transparent and unit-checked (evaluates to
# 5140.487 cm^-1, the standard value). A negative eigenvalue is an imaginary
# mode.
_HARTREE_J = 4.3597447222071e-18   # J
_BOHR_M = 0.529177210903e-10       # m
_AMU_KG = 1.66053906660e-27        # kg
_C_CM_S = 2.99792458e10            # cm/s
HESS_EIGVAL_TO_CM = (math.sqrt(_HARTREE_J / (_BOHR_M ** 2 * _AMU_KG))
                     / (2.0 * math.pi * _C_CM_S))

# Inner-sphere modes below this wavenumber (cm^-1) count as genuinely imaginary;
# the band between it and 0 is treated as numerical noise / a benign soft mode.
DEFAULT_IMAGINARY_THRESHOLD_CM = -50.0

# One perturb-and-reoptimise kick along the imaginary mode (plan Sec 5F) moves
# the mode's dominant atom by this much (Angstrom): enough to leave the saddle,
# small enough (a few percent of an M-L bond) that the re-optimisation stays in
# the nearby basin.
DEFAULT_PERTURB_STEP_A = 0.1


def read_xtb_hessian(path, n_atoms):
  """Parse an xtb ``.hessian`` file into a ``flex.double`` of shape (3N, 3N).

  Values are Hartree/Bohr^2, row-major, in the same atom order as the ``.xyz``
  fed to xtb.  Everything outside the ``$hessian ... $end`` numeric body
  (``$`` lines, blanks) is skipped.
  """
  vals = flex.double()
  with open(path) as fh:
    for line in fh:
      s = line.strip()
      if not s or s.startswith('$'):
        continue
      for tok in s.split():
        vals.append(float(tok))
  n = 3 * n_atoms
  if vals.size() != n * n:
    raise ValueError(f"hessian {path}: got {vals.size()} values, expected "
                     f"{n * n} for {n_atoms} atoms")
  vals.reshape(flex.grid(n, n))
  return vals


def read_orca_hessian(path, n_atoms):
  """Parse the ``$hessian`` block of an ORCA ``.hess`` file into a
  ``flex.double`` of shape (3N, 3N).

  A drop-in alternative to :func:`read_xtb_hessian`: same type, shape, units
  (Hartree/Bohr^2), row-major, in the input atom order.

  ORCA writes the (3N x 3N) Cartesian Hessian in column blocks.  After the
  ``$hessian`` line and a line with the matrix dimension M (= 3N) come the
  blocks: each opens with a header row of column indices (``0 1 2 3 4`` ...),
  then M rows whose leading integer is the row index followed by that block's
  values.  The block width is read from the header row, not assumed (typically
  5, sometimes 6).  A header row is distinguished from a data row by being all
  integers (data values are always in ``E`` notation).  The block ends at the
  next header row or ``$`` section (e.g. ``$vibrational_frequencies``).

  ORCA's raw Hessian is essentially symmetric; if asymmetry exceeds 1e-6 the
  result is symmetrised as ``0.5 (H + H^T)``, otherwise it is left as parsed.
  """
  n = 3 * n_atoms
  with open(path) as fh:
    lines = fh.readlines()
  i = 0
  while i < len(lines) and lines[i].strip() != '$hessian':
    i += 1
  if i >= len(lines):
    raise ValueError(f"hessian {path}: no $hessian block found")
  i += 1
  while i < len(lines) and not lines[i].strip():   # skip blanks before the dim
    i += 1
  if i >= len(lines):
    raise ValueError(f"hessian {path}: truncated $hessian block")
  dim = int(lines[i].split()[0])
  i += 1
  if dim != n:
    raise ValueError(f"hessian {path}: matrix dimension {dim}, expected "
                     f"{n} for {n_atoms} atoms")
  vals = flex.double(n * n, 0.0)
  seen = 0
  cols = None
  while i < len(lines):
    s = lines[i].strip()
    i += 1
    if not s:
      continue
    if s.startswith('$') or s.startswith('#'):     # next section: block done
      break
    toks = s.split()
    try:                                           # header row = all integers
      cols = [int(t) for t in toks]
      continue
    except ValueError:
      pass
    row = int(toks[0])                             # data row: index + values
    for k, t in enumerate(toks[1:]):
      vals[row * n + cols[k]] = float(t)
      seen += 1
  if seen != n * n:
    raise ValueError(f"hessian {path}: parsed {seen} values, expected "
                     f"{n * n} for {n_atoms} atoms")
  vals.reshape(flex.grid(n, n))
  asym = flex.max(flex.abs(vals - vals.matrix_transpose()))
  if asym > 1e-6:
    vals = (vals + vals.matrix_transpose()) * 0.5
  return vals


def _block(hessian, i, j):
  """The 3x3 Hessian sub-block for atoms (i, j) as ``scitbx.matrix.sqr``."""
  return matrix.sqr([hessian[(3 * i + a, 3 * j + b)]
                     for a in range(3) for b in range(3)])


def _projected_fc(block3, direction):
  """Seminario projection: ``sum_a lambda_a |u . v_a|`` over the symmetrised
  block's eigenpairs, where ``u`` is the unit *direction*.

  The off-diagonal Hessian block is not symmetric on its own
  (``H_ij = H_ji^T``); symmetrising avoids complex eigenvalues and makes the
  two per-atom estimates identical.  Returns a force constant in the block's
  units (Hartree/Bohr^2)."""
  sym = (block3 + block3.transpose()) * 0.5
  g = flex.double(sym.elems)
  g.reshape(flex.grid(3, 3))
  es = eigensystem.real_symmetric(g)
  values = es.values()
  vectors = es.vectors()
  u = direction.normalize()
  k = 0.0
  for a in range(3):
    v = matrix.col((vectors[3 * a], vectors[3 * a + 1], vectors[3 * a + 2]))
    k += values[a] * abs(u.dot(v))
  return k


def bond_force_constant(hessian, i, j, sites_cart):
  """Modified-Seminario M-L bond force constant in kcal/mol/Angstrom^2.

  *sites_cart* are the (optimised) Cartesian coordinates in Angstrom, in the
  Hessian's atom order."""
  u = matrix.col(sites_cart[j]) - matrix.col(sites_cart[i])
  k_hb2 = _projected_fc(-_block(hessian, i, j), u)
  return k_hb2 * HARTREE_BOHR2_TO_KCAL_A2


def angle_force_constant(hessian, a, b, c, sites_cart):
  """Modified-Seminario A-B-C angle force constant in kcal/mol/rad^2 (B is the
  central atom = the metal).  Returns ``None`` for a near-linear angle (no
  well-defined plane) or non-positive perpendicular curvature."""
  rb = matrix.col(sites_cart[b])
  vba = matrix.col(sites_cart[a]) - rb
  vbc = matrix.col(sites_cart[c]) - rb
  r_ba = vba.length()
  r_bc = vbc.length()
  if r_ba < 1e-6 or r_bc < 1e-6:
    return None
  u_ba = vba.normalize()
  u_bc = vbc.normalize()
  normal = u_ba.cross(u_bc)
  if normal.length() < 1e-4:          # near-linear: plane undefined
    return None
  normal = normal.normalize()
  u_perp_a = normal.cross(u_ba)       # perp to BA, in the ABC plane
  u_perp_c = u_bc.cross(normal)       # perp to BC, in the ABC plane
  k_ba = _projected_fc(-_block(hessian, b, a), u_perp_a)   # Hartree/Bohr^2
  k_bc = _projected_fc(-_block(hessian, b, c), u_perp_c)
  if k_ba <= 0.0 or k_bc <= 0.0:
    return None
  # Bond lengths in Bohr so k_perp * R^2 is an energy (Hartree); combine the
  # two perpendicular contributions (modified Seminario) -> Hartree/rad^2.
  r_ba_bohr = r_ba / BOHR_TO_ANG
  r_bc_bohr = r_bc / BOHR_TO_ANG
  inv = 1.0 / (k_ba * r_ba_bohr ** 2) + 1.0 / (k_bc * r_bc_bohr ** 2)
  return (1.0 / inv) * HARTREE_TO_KCAL


def bond_sigma(k_kcal_a2, temperature=DEFAULT_TEMPERATURE_K,
               sigma_min=DEFAULT_BOND_SIGMA_MIN,
               sigma_max=DEFAULT_BOND_SIGMA_MAX):
  """Bond restraint sigma (Angstrom) from a force constant via
  ``sigma = sqrt(k_B T / k)``, clamped to [sigma_min, sigma_max].
  Non-positive *k_kcal_a2* (numerical pathology) -> *sigma_max*."""
  if k_kcal_a2 is None or k_kcal_a2 <= 0.0:
    return sigma_max
  s = math.sqrt(KB_KCAL * temperature / k_kcal_a2)
  return min(max(s, sigma_min), sigma_max)


def angle_sigma(k_kcal_rad2, temperature=DEFAULT_TEMPERATURE_K,
                sigma_min=DEFAULT_ANGLE_SIGMA_MIN,
                sigma_max=DEFAULT_ANGLE_SIGMA_MAX):
  """Angle restraint sigma (degrees) from a force constant via
  ``sigma = sqrt(k_B T / k)``, clamped to [sigma_min, sigma_max].
  ``None``/non-positive *k_kcal_rad2* -> *sigma_max*."""
  if k_kcal_rad2 is None or k_kcal_rad2 <= 0.0:
    return sigma_max
  s_rad = math.sqrt(KB_KCAL * temperature / k_kcal_rad2)
  return min(max(math.degrees(s_rad), sigma_min), sigma_max)


def _classify_angle_sigmas(angle_ks, temperature,
                           sigma_max=DEFAULT_ANGLE_SIGMA_MAX):
  """Per-angle ``(sigma, origin)``, separating genuine curvatures from
  contaminated ones.

  The constrained-opt full Hessian is non-stationary, so its low-frequency
  donor-metal-donor bend projections are unreliable (some come back with
  negative curvature).  A curvature is trusted only when its force constant
  clears the clamp floor: positive and stiff enough to give a sigma below the
  ceiling.  Contaminated angles (``None`` / sub-floor k) take the site's LOOSEST
  trusted angle sigma (widen-only, ``seminario_fallback``); a site with no
  trusted angle drops its angles entirely.

  Returns ``(angle_sigmas {key: sigma_deg}, angle_origins {key: origin})``; a
  dropped angle is in neither dict.
  """
  k_floor = KB_KCAL * temperature / math.radians(sigma_max) ** 2
  trusted = {key: angle_sigma(k, temperature, sigma_max=sigma_max)
             for key, k in angle_ks.items()
             if k is not None and k > k_floor}
  fallback = max(trusted.values()) if trusted else None
  angle_sigmas, angle_origins = {}, {}
  for key in angle_ks:
    if key in trusted:
      angle_sigmas[key], angle_origins[key] = trusted[key], "hessian"
    elif fallback is not None:
      angle_sigmas[key], angle_origins[key] = fallback, "seminario_fallback"
  return angle_sigmas, angle_origins


def seminario_sigmas(hessian, metal_iseq, donor_iseqs, opt_sites_cart,
                     temperature=DEFAULT_TEMPERATURE_K, log=None):
  """Modified-Seminario sigmas for every M-L bond and L-M-L angle of a site.

  Parameters
  ----------
  hessian : flex.double (3N, 3N)
      From :func:`read_xtb_hessian`; Hartree/Bohr^2, at *opt_sites_cart*.
  metal_iseq : int
  donor_iseqs : sequence of int
  opt_sites_cart : scitbx.array_family.flex.vec3_double
      Optimised geometry (Angstrom), Hessian atom order.
  temperature : float
  log : file-like or None

  Returns
  -------
  libtbx.group_args
      ``bond_sigmas`` {donor_iseq: sigma_A}, ``angle_sigmas``
      {frozenset({i, j}): sigma_deg}, ``angle_origins`` {frozenset({i, j}):
      ``hessian`` | ``seminario_fallback``} parallel to *angle_sigmas*, plus the
      raw force constants ``bond_ks`` / ``angle_ks`` (for provenance /
      reporting).
  """
  if log is None:
    log = null_out()
  bond_sigmas, bond_ks = {}, {}
  for d in donor_iseqs:
    k = bond_force_constant(hessian, metal_iseq, d, opt_sites_cart)
    bond_ks[d] = k
    bond_sigmas[d] = bond_sigma(k, temperature)

  angle_ks = {}
  donors = list(donor_iseqs)
  for x in range(len(donors)):
    for y in range(x + 1, len(donors)):
      i, j = donors[x], donors[y]
      angle_ks[frozenset((i, j))] = angle_force_constant(
        hessian, i, metal_iseq, j, opt_sites_cart)
  angle_sigmas, angle_origins = _classify_angle_sigmas(angle_ks, temperature)

  print(f"    Seminario sigmas: {len(bond_sigmas)} bond, "
        f"{len(angle_sigmas)} angle emitted of {len(angle_ks)} "
        f"(T={temperature:.1f} K)", file=log)
  return group_args(bond_sigmas=bond_sigmas, angle_sigmas=angle_sigmas,
                    angle_origins=angle_origins, bond_ks=bond_ks,
                    angle_ks=angle_ks)


def _mass_weighted_inner_block(hessian, iseqs, masses):
  """The mass-weighted metal+donor Hessian sub-block and its per-DOF
  1/sqrt(mass).

  Gathers the 3m x 3m Cartesian block for atoms *iseqs* (m = len(iseqs)) from
  the full *hessian* by a row-major flex ``select``, symmetrises it, and
  mass-weights (``D H D`` with ``D = diag(1/sqrt(m))``, an outer-product
  elementwise multiply). Returns
  ``(hw [flex dim x dim], inv_sqrt_m [flex dim])``.

  Only the metal+donor block is used (plan Sec 5F): those atoms were free during
  the constrained opt, so this block is the exact second-order condition for
  that constrained minimum. The clamped scaffold's non-stationarity lives only
  in the discarded cross/scaffold blocks and cannot enter here; and because the
  scaffold is clamped, the inner cluster's rigid translation/rotation is
  resisted, so there are no trivial zero modes to project out."""
  n3 = hessian.all()[0]
  dim = 3 * len(iseqs)
  dof = [3 * int(ip) + a for ip in iseqs for a in range(3)]     # 3m DOF indices
  flat = flex.size_t([r * n3 + c for r in dof for c in dof])  # row-major picks
  sub = hessian.as_1d().select(flat)
  sub.reshape(flex.grid(dim, dim))
  sub = (sub + sub.matrix_transpose()) * 0.5                    # symmetrise
  inv_sqrt_m = flex.double([1.0 / math.sqrt(masses[int(ip)])
                            for ip in iseqs for _ in range(3)])
  wcol = inv_sqrt_m.deep_copy(); wcol.reshape(flex.grid(dim, 1))
  wrow = inv_sqrt_m.deep_copy(); wrow.reshape(flex.grid(1, dim))
  return sub * wcol.matrix_multiply(wrow), inv_sqrt_m


def _inner_eigensystem(hessian, iseqs, masses):
  """Eigenpairs of the mass-weighted metal+donor Hessian sub-block, plus the
  per-DOF 1/sqrt(mass). Diagonalises once so the imaginary-mode check and the
  perturbation can share it (both need the same eigenpairs).

  Returns ``group_args(vals, vecs, inv_sqrt_m, dim)``: *vals* the eigenvalues
  (omega^2), *vecs* the flat eigenvector matrix, *inv_sqrt_m* the 1/sqrt(mass)
  per DOF, *dim* = 3m."""
  hw, inv_sqrt_m = _mass_weighted_inner_block(hessian, iseqs, masses)
  es = eigensystem.real_symmetric(hw)
  return group_args(vals=es.values(), vecs=es.vectors(),
                    inv_sqrt_m=inv_sqrt_m, dim=3 * len(iseqs))


def _frequencies_cm(vals):
  """Sorted harmonic wavenumbers (cm^-1) from eigenvalues (omega^2); a NEGATIVE
  value marks an imaginary mode (omega^2 < 0)."""
  cm = flex.sqrt(flex.abs(vals)) * HESS_EIGVAL_TO_CM            # |nu| in cm^-1
  neg = vals < 0.0
  cm.set_selected(neg, -cm.select(neg))  # imaginary -> < 0
  return sorted(cm)


def inner_sphere_frequencies(hessian, iseqs, masses):
  """Harmonic wavenumbers (cm^-1) of the metal+donor Hessian sub-block, sorted
  ascending; a NEGATIVE value is an imaginary mode (omega^2 < 0).

  Diagonalises the mass-weighted sub-block (:func:`_inner_eigensystem`).
  Mass-weighting does not change WHICH modes are imaginary (Sylvester's law of
  inertia); it sets the frequency scale so the threshold is physical and a real
  imaginary mode is reportable in cm^-1."""
  return _frequencies_cm(_inner_eigensystem(hessian, iseqs, masses).vals)


def imaginary_mode_perturbation(hessian, metal_iseq, donor_iseqs, masses,
                                sites_cart, step=DEFAULT_PERTURB_STEP_A,
                                eigen=None):
  """Perturb *sites_cart* along the most-imaginary inner-sphere mode (plan Sec
  5F).

  Returns a new ``flex.vec3_double`` where only the metal+donor atoms move,
  along the imaginary mode's Cartesian displacement (its mass-weighted
  eigenvector, un-weighted per DOF), scaled so the mode's DOMINANT atom moves
  *step* Angstrom. Deterministic (the eigenvector is fixed by the Hessian), so
  runs stay reproducible; the sign is arbitrary (either direction escapes the
  saddle).

  *eigen*, when supplied (the ``eigen`` field of an :func:`imaginary_mode_check`
  result), reuses that diagonalisation instead of repeating it."""
  iseqs = [metal_iseq] + list(donor_iseqs)
  if eigen is None:
    eigen = _inner_eigensystem(hessian, iseqs, masses)
  vals, vecs, inv_sqrt_m, dim = (eigen.vals, eigen.vecs,
                                 eigen.inv_sqrt_m, eigen.dim)
  i_min = min(range(dim), key=lambda i: vals[i])               # most-imaginary
  q = flex.double([vecs[i_min * dim + k] for k in range(dim)])
  dx = q * inv_sqrt_m  # mass-weighted -> Cartesian
  disp = [matrix.col((dx[3 * p], dx[3 * p + 1], dx[3 * p + 2]))
          for p in range(len(iseqs))]
  d_max = max(v.length() for v in disp)
  out = sites_cart.deep_copy()
  if d_max < 1e-12:
    return out
  scale = step / d_max  # dominant atom -> step A
  for ip, v in zip(iseqs, disp):
    out[ip] = tuple(matrix.col(out[ip]) + v * scale)
  return out


def imaginary_mode_check(hessian, metal_iseq, donor_iseqs, masses,
                         threshold_cm=DEFAULT_IMAGINARY_THRESHOLD_CM, log=None):
  """Whether the inner sphere (metal + donors) sits at a true local minimum.

  Diagonalises the mass-weighted metal+donor Hessian sub-block
  (:func:`_inner_eigensystem`) and treats any mode below *threshold_cm* as a
  genuine imaginary mode (the band above it is numerical noise / a benign soft
  mode). Returns group_args(stable, frequencies, imaginary, eigen): *stable* is
  True when no mode is below threshold, *frequencies* the full sorted list
  (cm^-1), *imaginary* the offending ones, *eigen* the diagonalisation for
  :func:`imaginary_mode_perturbation` to reuse. The trigger for the
  perturb-and-reoptimise step (plan Sec 5F)."""
  if log is None:
    log = null_out()
  eigen = _inner_eigensystem(hessian, [metal_iseq] + list(donor_iseqs), masses)
  freqs = _frequencies_cm(eigen.vals)
  imaginary = [f for f in freqs if f < threshold_cm]
  if imaginary:
    modes = ", ".join(f"{-f:.0f}i" for f in imaginary)
    print(f"    inner-sphere imaginary mode(s) at {modes} cm^-1 "
          f"(threshold {threshold_cm:.0f})", file=log)
  return group_args(stable=not imaginary, frequencies=freqs,
                    imaginary=imaginary, eigen=eigen)
