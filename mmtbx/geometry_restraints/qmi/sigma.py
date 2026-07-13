"""Modified-Seminario force constants -> restraint sigmas (plan Sec 5F).

Reads the Cartesian Hessian xTB writes (``--hess``, a ``<namespace>.hessian``
file, Hartree/Bohr^2, full 3N x 3N) at the optimised geometry and derives a
per-restraint force constant by the modified-Seminario method (Allen, Payne &
Cole 2018, improving Seminario 1996):

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

from __future__ import absolute_import, division, print_function

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
DEFAULT_TEMPERATURE_K = 298.15

# Sigma clamps (plan: sigmas are floors to widen; never emit a false-tight
# restraint from numerical noise).
DEFAULT_BOND_SIGMA_MIN = 0.005    # A
DEFAULT_BOND_SIGMA_MAX = 0.20     # A
DEFAULT_ANGLE_SIGMA_MIN = 0.5     # deg
DEFAULT_ANGLE_SIGMA_MAX = 15.0    # deg


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
      {frozenset({i, j}): sigma_deg}, plus the raw force constants
      ``bond_ks`` / ``angle_ks`` (for provenance / reporting).
  """
  if log is None:
    log = null_out()
  bond_sigmas, bond_ks = {}, {}
  for d in donor_iseqs:
    k = bond_force_constant(hessian, metal_iseq, d, opt_sites_cart)
    bond_ks[d] = k
    bond_sigmas[d] = bond_sigma(k, temperature)

  angle_sigmas, angle_ks = {}, {}
  donors = list(donor_iseqs)
  for x in range(len(donors)):
    for y in range(x + 1, len(donors)):
      i, j = donors[x], donors[y]
      k = angle_force_constant(hessian, i, metal_iseq, j, opt_sites_cart)
      key = frozenset((i, j))
      angle_ks[key] = k
      angle_sigmas[key] = angle_sigma(k, temperature)

  print(f"    Seminario sigmas: {len(bond_sigmas)} bond, "
        f"{len(angle_sigmas)} angle (T={temperature:.1f} K)", file=log)
  return group_args(bond_sigmas=bond_sigmas, angle_sigmas=angle_sigmas,
                    bond_ks=bond_ks, angle_ks=angle_ks)
