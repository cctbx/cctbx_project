"""Regression tests for mmtbx.geometry_restraints.qmi (QM metal restraints).

Covers the qmi-specific, deterministic logic: the sigma math (Hessian parse +
modified Seminario + sigma conversion), first-shell detection, restraint
assembly + .edits emission, the frozen/free partition, and the metals helpers.

Run: mmtbx.python mmtbx/regression/tst_qmi.py
"""

from __future__ import absolute_import, division, print_function

import math
import os
import tempfile
from io import StringIO

import iotbx.pdb
import libtbx.phil
import mmtbx.model
from iotbx.data_manager import DataManager
from libtbx import group_args
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out, Sorry
from scitbx import matrix
from scitbx.array_family import flex

from mmtbx.programs.qmi import Program, master_phil_str
from mmtbx.geometry_restraints import orca_manager, xtb_manager
from mmtbx.geometry_restraints.qmi import (
  qm, restraints, sigma, metals, protonation, shape, bvs, sanity, provenance,
  util)
from mmtbx.geometry_restraints.qmi import electronic as el
from mmtbx.geometry_restraints.qmi import donors as D
from mmtbx.geometry_restraints.qmi import metal_protonation as mp
from mmtbx.geometry_restraints.qmi.builder import QMRestraintsBuilder


# A minimal Zn(II) N3O site: Zn + 3 His N donors + a catalytic water O, plus a
# His ring C (~2.8 A, must be excluded by the element filter) and an H (must be
# excluded).  Geometry is nominal but exact: His94 NE2 on +x at 2.0 A, His96
# NE2 on +y at 2.05 A -> a 90 deg N-Zn-N angle to check against.
_ZN_SITE_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  NE2 HIS A  94       2.000   0.000   0.000  1.00 10.00           N
ATOM      3  NE2 HIS A  96       0.000   2.050   0.000  1.00 10.00           N
ATOM      4  ND1 HIS A 119       0.000   0.000   2.050  1.00 10.00           N
HETATM    5  O   HOH A 263      -1.400  -1.400   0.000  1.00 10.00           O
ATOM      6  CD2 HIS A  94       2.500   0.900   0.900  1.00 10.00           C
ATOM      7  HE2 HIS A  94       1.900   0.500   0.000  1.00 10.00           H
"""

# Metal + a near N donor (free) + a far water O (> free_radius, frozen) + a near
# H (free), for the frozen/free partition test.
_FREEZE_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A   1       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  NE2 HIS A   2       2.000   0.000   0.000  1.00 10.00           N
HETATM    3  O   HOH A   3       0.000   0.000  10.000  1.00 10.00           O
ATOM      4  HE2 HIS A   2       1.500   0.000   0.000  1.00 10.00           H
"""

# Metal + a coordinating water O (~2.1 A), a non-coordinating water O in the
# relax-zone annulus (~4 A, within free_radius), and a non-coordinating water O
# beyond free_radius (~7 A), for the freeze_waters mode test.
_WATER_MODE_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A   1       0.000   0.000   0.000  1.00 10.00          ZN
HETATM    2  O   HOH A   2       2.100   0.000   0.000  1.00 10.00           O
HETATM    3  O   HOH A   3       4.000   0.000   0.000  1.00 10.00           O
HETATM    4  O   HOH A   4       7.000   0.000   0.000  1.00 10.00           O
"""


# A Zn + an unprotonated donor (no H) for the validate protonation guard.
_NOH_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  NE2 HIS A  94       2.000   0.000   0.000  1.00 10.00           N
"""

# Zn + a HIP His (both ring N protonated) coordinating via NE2 + a coordinating
# Cys, for the metal-donor deprotonation test.
_HIS_HIP_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  NE2 HIS A  94       2.000   0.000   0.000  1.00 10.00           N
ATOM      3  CE1 HIS A  94       2.600   1.100   0.000  1.00 10.00           C
ATOM      4  CD2 HIS A  94       2.600  -1.100   0.000  1.00 10.00           C
ATOM      5  ND1 HIS A  94       3.900   0.700   0.000  1.00 10.00           N
ATOM      6  CG  HIS A  94       3.900  -0.700   0.000  1.00 10.00           C
ATOM      7  HE2 HIS A  94       1.300   0.700   0.000  1.00 10.00           H
ATOM      8  HD1 HIS A  94       4.300   1.600   0.000  1.00 10.00           H
ATOM      9  SG  CYS A  99       0.000   2.300   0.000  1.00 10.00           S
ATOM     10  HG  CYS A  99       0.500   3.000   0.000  1.00 10.00           H
"""

# Zn + a wrong-tautomer His: the proton sits only on the coordinating N (NE2),
# so it must be relocated to the partner N (ND1).
_HIS_TAUT_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  NE2 HIS A  94       2.000   0.000   0.000  1.00 10.00           N
ATOM      3  CE1 HIS A  94       2.600   1.100   0.000  1.00 10.00           C
ATOM      4  CD2 HIS A  94       2.600  -1.100   0.000  1.00 10.00           C
ATOM      5  ND1 HIS A  94       3.900   0.700   0.000  1.00 10.00           N
ATOM      6  CG  HIS A  94       3.900  -0.700   0.000  1.00 10.00           C
ATOM      7  HE2 HIS A  94       1.300   0.700   0.000  1.00 10.00           H
"""

# His with a proton but no metal: deprotonation must be a no-op.
_NOMETAL_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
ATOM      1  NE2 HIS A  94       2.000   0.000   0.000  1.00 10.00           N
ATOM      2  HE2 HIS A  94       1.300   0.700   0.000  1.00 10.00           H
"""

# Zn + a bidentate coordinating Asp carboxylate (both OD within cutoff), OD2
# protonated. Exercises the anionic strip AND the bidentate dedup (both O keys
# map to the same HD2, so it must be removed exactly once).
_ASP_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  OD1 ASP A  50       2.000   0.200   0.000  1.00 10.00           O
ATOM      3  OD2 ASP A  50       2.000  -0.200   0.000  1.00 10.00           O
ATOM      4  CG  ASP A  50       3.000   0.000   0.000  1.00 10.00           C
ATOM      5  HD2 ASP A  50       2.400   1.000   0.000  1.00 10.00           H
"""

# Zn + a coordinating Cys whose thiol proton is a DEUTERIUM (DG); the strip rule
# must remove DG as well as HG.
_CYS_D_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  SG  CYS A  99       0.000   2.300   0.000  1.00 10.00           S
ATOM      3  DG  CYS A  99       0.500   3.000   0.000  1.00 10.00           D
"""

# Mg (a hard ion) coordinating a Cys and an Asp: the metal-induced Cys thiolate
# must NOT form (HG stays), but the pH-driven Asp carboxylate must (HD2 goes).
_MG_DONORS_PDB = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1 MG    MG A 262       0.000   0.000   0.000  1.00 10.00          MG
ATOM      2  SG  CYS A  99       2.300   0.000   0.000  1.00 10.00           S
ATOM      3  HG  CYS A  99       2.900   0.600   0.000  1.00 10.00           H
ATOM      4  OD1 ASP A  50       0.000   2.000   0.000  1.00 10.00           O
ATOM      5  OD2 ASP A  50       0.000   2.400   0.300  1.00 10.00           O
ATOM      6  CG  ASP A  50       0.000   3.000   0.000  1.00 10.00           C
ATOM      7  HD2 ASP A  50       0.500   2.500   0.000  1.00 10.00           H
"""

# Zn + a bidentate Asp (both OD within cutoff, ~60 deg bite) + two His N, for
# the carboxylate-regularisation shape test.
_CARBOX_SHELL_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  OD1 ASP A  50       1.905   1.100   0.000  1.00 10.00           O
ATOM      3  OD2 ASP A  50       1.905  -1.100   0.000  1.00 10.00           O
ATOM      4  CG  ASP A  50       2.800   0.000   0.000  1.00 10.00           C
ATOM      5  NE2 HIS A  94       0.000   2.000   0.000  1.00 10.00           N
ATOM      6  NE2 HIS A  96       0.000   0.000   2.000  1.00 10.00           N
"""

# Zn + an anisobidentate Asp (OD1 at 2.1 A, OD2 at 2.7 A: monodentate + a long
# contact, not a chelate) + two His N.
_CARBOX_ANISO_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  OD1 ASP A  50       1.819   1.050   0.000  1.00 10.00           O
ATOM      3  OD2 ASP A  50       2.338  -1.350   0.000  1.00 10.00           O
ATOM      4  CG  ASP A  50       2.800   0.000   0.000  1.00 10.00           C
ATOM      5  NE2 HIS A  94       0.000   2.000   0.000  1.00 10.00           N
ATOM      6  NE2 HIS A  96       0.000   0.000   2.000  1.00 10.00           N
"""

# Zn + a C-terminal Asp chelating via its BACKBONE carboxylate (O/OXT, symmetric)
# + two His N. The coordinating pair is the main-chain O/OXT, not the sidechain.
_CTERM_CARBOX_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  O   ASP A  50       1.905   1.100   0.000  1.00 10.00           O
ATOM      3  OXT ASP A  50       1.905  -1.100   0.000  1.00 10.00           O
ATOM      4  C   ASP A  50       2.800   0.000   0.000  1.00 10.00           C
ATOM      5  NE2 HIS A  94       0.000   2.000   0.000  1.00 10.00           N
ATOM      6  NE2 HIS A  96       0.000   0.000   2.000  1.00 10.00           N
"""

# Zn + a bare coordinating water, for the cation-aware water-H placement.
_WATER_METAL_PDB = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
HETATM    2  O   HOH A 300       2.100   0.000   0.000  1.00 10.00           O
"""

# Cu + a coordinating backbone amide (ALA N-H, ~2 A) + a coordinating N-terminal
# amine (GLY -NH3+, ~2 A), for the ATCUN deprotonation (amidate + -NH2 amine).
_ATCUN_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 CU    CU A 262       0.000   0.000   0.000  1.00 10.00          CU
ATOM      2  N   ALA A   3       2.000   0.000   0.000  1.00 10.00           N
ATOM      3  H   ALA A   3       2.600   0.700   0.000  1.00 10.00           H
ATOM      4  N   GLY A   1       0.000   2.000   0.000  1.00 10.00           N
ATOM      5  H1  GLY A   1       0.600   2.600   0.000  1.00 10.00           H
ATOM      6  H2  GLY A   1      -0.600   2.600   0.000  1.00 10.00           H
ATOM      7  H3  GLY A   1       0.000   2.600   0.600  1.00 10.00           H
"""

# The same backbone amide near a Zn (not an amidate-inducing ion): the amide H
# must be kept (metal-gating test).
_ZN_AMIDE_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  N   ALA A   3       2.000   0.000   0.000  1.00 10.00           N
ATOM      3  H   ALA A   3       2.600   0.700   0.000  1.00 10.00           H
"""

# Zn + a MONODENTATE Asp: OD1 coordinates (2.0 A), OD2 dangles beyond the cutoff
# (3.8 A), for the far-O freeze.
_CARBOX_MONO_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  OD1 ASP A  50       2.000   0.000   0.000  1.00 10.00           O
ATOM      3  OD2 ASP A  50       3.500   1.500   0.000  1.00 10.00           O
ATOM      4  CG  ASP A  50       3.000   0.000   0.000  1.00 10.00           C
ATOM      5  NE2 HIS A  94       0.000   2.000   0.000  1.00 10.00           N
"""

# Zn + a His N donor (2 A) + a non-donor water O (4 A), for the coordination-mode
# check: the test moves the metal-relative geometry to make one leave, one enter.
_MODE_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
ATOM      2  NE2 HIS A  94       2.000   0.000   0.000  1.00 10.00           N
HETATM    3  O   HOH A 300       4.000   0.000   0.000  1.00 10.00           O
"""


def _model(pdb_str):
  """mmtbx.model.manager from an inline PDB string (no processing)."""
  pdb_in = iotbx.pdb.input(source_info=None, lines=pdb_str.splitlines())
  return mmtbx.model.manager(model_input=pdb_in)


# ---------------------------------------------------------------------------
# metals
# ---------------------------------------------------------------------------
def exercise_metals():
  assert metals.is_metal("ZN")
  assert metals.is_metal("Zn")
  assert metals.is_metal(" zn ")
  assert metals.is_metal("Mn")
  assert not metals.is_metal("C")
  assert not metals.is_metal("O")
  model = _model(_ZN_SITE_PDB)
  assert metals.count_metals(model) == 1
  found = metals.metal_atoms(model)
  assert len(found) == 1
  assert found[0].element.strip().capitalize() == "Zn"
  assert len(metals.metal_atoms(model, metals={"Fe"})) == 0   # element filter


# ---------------------------------------------------------------------------
# sigma: Hessian parse
# ---------------------------------------------------------------------------
def exercise_hessian_parse():
  # 2-atom system -> 6x6 = 36 values, row-major.
  vals = [float(i) for i in range(36)]
  tmp = tempfile.mkdtemp()
  path = os.path.join(tmp, "probe.hessian")
  with open(path, "w") as fh:
    fh.write("$hessian\n")
    for i in range(0, 36, 5):
      fh.write("  " + "  ".join(f"{v:.10f}" for v in vals[i:i + 5]) + "\n")
    fh.write("$end\n")
  h = sigma.read_xtb_hessian(path, 2)
  assert h.focus() == (6, 6), h.focus()
  assert approx_equal(h[(0, 0)], 0.0)
  assert approx_equal(h[(0, 5)], 5.0)
  assert approx_equal(h[(5, 5)], 35.0)
  try:
    sigma.read_xtb_hessian(path, 3)     # wrong atom count must raise
  except ValueError:
    pass
  else:
    raise AssertionError("expected ValueError on size mismatch")


# ---------------------------------------------------------------------------
# sigma: modified-Seminario bond force constant (analytic recovery)
# ---------------------------------------------------------------------------
def exercise_seminario_bond():
  # A single stretch of constant k along +x between atoms 0 and 1 gives the
  # Hessian off-diagonal block H_01 = -k (u u^T), u = x_hat, i.e. only
  # H[0,3] = -k.  The Seminario projection onto u must recover k (Hartree/
  # Bohr^2), returned in kcal/mol/A^2.
  k_hb2 = 0.35
  vals = flex.double(36, 0.0)
  vals[0 * 6 + 3] = -k_hb2             # H[0,3]
  vals[3 * 6 + 0] = -k_hb2             # H[3,0] (keep the full matrix symmetric)
  vals.reshape(flex.grid(6, 6))
  sites = flex.vec3_double([(0.0, 0.0, 0.0), (2.0, 0.0, 0.0)])
  k = sigma.bond_force_constant(vals, 0, 1, sites)
  assert approx_equal(k, k_hb2 * sigma.HARTREE_BOHR2_TO_KCAL_A2, eps=1e-4), k
  # A bond along +y sees no curvature in this (x-only) block -> ~0.
  sites_y = flex.vec3_double([(0.0, 0.0, 0.0), (0.0, 2.0, 0.0)])
  assert approx_equal(
    sigma.bond_force_constant(vals, 0, 1, sites_y), 0.0, eps=1e-6)


def exercise_seminario_angle_guards():
  h = flex.double(81, 0.0)             # 3 atoms -> 9x9, all zero
  h.reshape(flex.grid(9, 9))
  # near-linear a-b-c (b central) -> plane undefined -> None
  collinear = flex.vec3_double([(2.0, 0.0, 0.0), (0.0, 0.0, 0.0),
                                (-2.0, 0.0, 0.0)])
  assert sigma.angle_force_constant(h, 0, 1, 2, collinear) is None
  # bent but zero perpendicular curvature -> None
  bent = flex.vec3_double([(2.0, 0.0, 0.0), (0.0, 0.0, 0.0), (0.0, 2.0, 0.0)])
  assert sigma.angle_force_constant(h, 0, 1, 2, bent) is None


# ---------------------------------------------------------------------------
# sigma: force-constant -> sigma conversion + clamps
# ---------------------------------------------------------------------------
def exercise_sigma_conversion():
  k = 100.0
  expect = math.sqrt(sigma.KB_KCAL * sigma.DEFAULT_TEMPERATURE_K / k)
  s = sigma.bond_sigma(k, sigma_min=0.0, sigma_max=1.0)
  assert approx_equal(s, expect, eps=1e-6), (s, expect)
  # temperature dependence: sigma proportional to sqrt(T)
  s2 = sigma.bond_sigma(k, temperature=4 * sigma.DEFAULT_TEMPERATURE_K,
                        sigma_min=0.0, sigma_max=1.0)
  assert approx_equal(s2, 2.0 * expect, eps=1e-6), (s2, expect)
  # non-positive / None curvature -> loose cap
  assert sigma.bond_sigma(0.0, sigma_max=0.2) == 0.2
  assert sigma.bond_sigma(-5.0, sigma_max=0.2) == 0.2
  assert sigma.angle_sigma(None, sigma_max=15.0) == 15.0
  # clamps
  assert sigma.bond_sigma(1e9, sigma_min=0.01, sigma_max=0.2) == 0.01
  assert sigma.bond_sigma(1e-9, sigma_min=0.01, sigma_max=0.2) == 0.2


def exercise_classify_angle_sigmas():
  """Contaminated angles (None / sub-floor k, from the non-stationary Hessian)
  take the site's LOOSEST trusted angle sigma (widen-only, seminario_fallback);
  a site with no trusted angle drops its angles; trusted angles keep their own
  hessian sigma."""
  T = sigma.DEFAULT_TEMPERATURE_K
  ka, kb, kc, kd = (frozenset((0, 1)), frozenset((0, 2)),
                    frozenset((1, 2)), frozenset((0, 3)))
  # a stiff (small sigma), b soft-but-trusted (larger sigma), c None, d sub-floor
  s, o = sigma._classify_angle_sigmas({ka: 100.0, kb: 10.0, kc: None, kd: 1.0}, T)
  assert o[ka] == "hessian" and approx_equal(s[ka], sigma.angle_sigma(100.0, T))
  assert o[kb] == "hessian" and approx_equal(s[kb], sigma.angle_sigma(10.0, T))
  # contaminated c/d widen to the site's loosest trusted angle (b), never tighten
  assert o[kc] == "seminario_fallback" and s[kc] == s[kb]
  assert o[kd] == "seminario_fallback" and s[kd] == s[kb]
  assert s[kb] == max(s[ka], s[kb])
  # no trusted angle -> every angle dropped (present in neither dict)
  s2, o2 = sigma._classify_angle_sigmas({ka: None, kb: -3.0}, T)
  assert s2 == {} and o2 == {}


def exercise_emit_angles():
  """Angle-restraint gate (builder._emit_angles): 'auto' emits when the emitted
  M-L bonds fail to pin the metal in 3D -- <=2 bonds (metal on a circle/line), or
  >=3 bonds whose directions are (near-)coplanar (a planar shell with a soft
  out-of-plane mode). A 3D-spanning shell needs no angle. 'always'/'never'
  override."""
  def dirs(*vs):
    return [matrix.col(v).normalize() for v in vs]
  auto = QMRestraintsBuilder(angle_restraints="auto", logger=null_out())
  always = QMRestraintsBuilder(angle_restraints="always", logger=null_out())
  never = QMRestraintsBuilder(angle_restraints="never", logger=null_out())

  linear = dirs((1, 0, 0), (-1, 0, 0))
  trigonal = dirs((1, 0, 0), (-0.5, 0.866, 0), (-0.5, -0.866, 0))
  square_planar = dirs((1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0))
  tetrahedral = dirs((1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1))
  octahedral = dirs((1, 0, 0), (-1, 0, 0), (0, 1, 0),
                    (0, -1, 0), (0, 0, 1), (0, 0, -1))

  # auto: coplanar / under-determined -> emit; a 3D-spanning shell -> no angle.
  assert auto._emit_angles(linear) is True             # <=2: metal on a circle
  assert auto._emit_angles(trigonal) is True           # 3 coplanar bonds
  assert auto._emit_angles(square_planar) is True        # 4 coplanar bonds
  assert auto._emit_angles(tetrahedral) is False         # spans 3D
  assert auto._emit_angles(octahedral) is False          # spans 3D
  for ds in (trigonal, tetrahedral):                   # modes ignore geometry
    assert always._emit_angles(ds) is True
    assert never._emit_angles(ds) is False


def exercise_coordination_number():
  """builder._coordination_number: the chemical (denticity) count fed to the
  predictor is effective cn plus one per symmetric bidentate carboxylate (an
  anisobidentate far O stays out, already excluded from cn)."""
  b = QMRestraintsBuilder(logger=null_out())
  # 2 His + monodentate Asp: no chelate -> denticity == effective cn.
  assert b._coordination_number(group_args(cn=3, n_bidentate_carboxylate=0)) == 3
  # 2 His + chelating Glu: cn=3 (2 N + bisector), one chelate -> denticity 4.
  assert b._coordination_number(group_args(cn=3, n_bidentate_carboxylate=1)) == 4
  # two chelates -> add both back.
  assert b._coordination_number(group_args(cn=4, n_bidentate_carboxylate=2)) == 6


def exercise_imaginary_mode():
  # 3 atoms: metal Fe (0), donor N (1), scaffold C (2). Diagonal Hessian; only
  # the metal+donor sub-block is diagonalised, so the scaffold's huge negative
  # curvature must NOT appear.
  n = 3
  H = flex.double(flex.grid(3 * n, 3 * n), 0.0)
  diag = [0.5, 0.5, 0.5,        # Fe (0): stable
          0.5, 0.5, -0.02,      # N  (1): one imaginary direction
          -10., -10., -10.]     # C  (2): unstable, but NOT inner -> excluded
  for i in range(3 * n):
    H[(i, i)] = diag[i]
  masses = [55.847, 14.007, 12.011]                # Fe, N, C (amu)
  chk = sigma.imaginary_mode_check(H, 0, [1], masses, log=null_out())
  assert not chk.stable and len(chk.imaginary) == 1
  # magnitude: -0.02 Hartree/Bohr^2 on 14.007 amu -> ~ -194 cm^-1
  expect = -sigma.HESS_EIGVAL_TO_CM * math.sqrt(0.02 / 14.007)
  assert approx_equal(chk.imaginary[0], expect, eps=1.0), (chk.imaginary[0], expect)
  # scaffold's -10 (which would be ~ -4700 cm^-1) is excluded by the sub-block
  assert min(chk.frequencies) > -300.0, chk.frequencies
  # sign convention: mass-weighting preserves the count of imaginary modes
  assert sum(1 for f in chk.frequencies if f < 0) == 1
  # perturb along the imaginary mode: only its dominant atom (N, along z) moves,
  # by `step` A; the metal (0, inner but off the mode) and the scaffold (2, not
  # inner) are untouched. Deterministic.
  sites = flex.vec3_double([(0, 0, 0), (2.0, 0, 0), (0, 2.0, 0)])
  out = sigma.imaginary_mode_perturbation(H, 0, [1], masses, sites, step=0.1)
  assert approx_equal(out[1], (2.0, 0.0, 0.1), eps=1e-6), out[1]
  assert approx_equal(out[0], sites[0], eps=1e-9)
  assert approx_equal(out[2], sites[2], eps=1e-9)
  # a fully stable inner sphere -> no imaginary mode
  H2 = flex.double(flex.grid(3 * n, 3 * n), 0.0)
  for i in range(3 * n):
    H2[(i, i)] = 0.5
  chk2 = sigma.imaginary_mode_check(H2, 0, [1], masses, log=null_out())
  assert chk2.stable and not chk2.imaginary


def exercise_perturb_reopt():
  """Perturb-and-reoptimise orchestration (plan Sec 5F) with the QM layer stubbed
  (no xTB): stable -> hessian; imaginary-then-stable -> reopt then hessian;
  imaginary past the attempt budget -> conservative_fallback. Multiple attempts
  escalate the kick and use a distinct namespace each time."""
  b1 = QMRestraintsBuilder(perturb_attempts=1, logger=null_out())
  b3 = QMRestraintsBuilder(perturb_attempts=3, logger=null_out())
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  hpre, opre, steps = [], [], []
  calls = {"reopt": 0}
  stability = iter(())
  saved = (qm.run_hessian, qm.run_constrained_opt, sigma.read_xtb_hessian,
           sigma.seminario_sigmas, sigma.imaginary_mode_check,
           sigma.imaginary_mode_perturbation)
  try:
    def fake_hess(*a, **k):
      hpre.append(k.get("preamble")); return "hess_path"
    def fake_opt(*a, **k):
      calls["reopt"] += 1; opre.append(k.get("preamble")); return sites, None
    def fake_perturb(*a, **k):
      steps.append(k.get("step")); return sites
    qm.run_hessian = fake_hess
    qm.run_constrained_opt = fake_opt
    sigma.imaginary_mode_perturbation = fake_perturb
    sigma.read_xtb_hessian = lambda *a, **k: None
    sigma.seminario_sigmas = lambda *a, **k: group_args(
      bond_sigmas={1: 0.02}, angle_sigmas={}, angle_origins={})
    sigma.imaginary_mode_check = lambda *a, **k: group_args(
      stable=next(stability), frequencies=[], imaginary=[], eigen=None)

    def run(builder):
      hpre[:] = []; opre[:] = []; steps[:] = []; calls["reopt"] = 0
      return builder._sigmas_with_reopt(1, model, sites, 0, [1, 2, 3, 4], [],
                                        0, 1, None)
    stability = iter([True])                               # stable first
    _s, bs, _as, _ao, origin, npert = run(b1)
    assert origin == "hessian" and bs == {1: 0.02}
    assert calls["reopt"] == 0 and npert == 0
    stability = iter([False, True])                        # imaginary -> reopt -> ok
    _s, bs, _as, _ao, origin, npert = run(b1)
    assert origin == "hessian" and calls["reopt"] == 1 and npert == 1
    assert hpre == ["qmi_hess", "qmi_reopt1_hess"], hpre   # distinct namespaces
    stability = iter([False, False])                       # past the 1-attempt budget
    _s, bs, _as, _ao, origin, npert = run(b1)
    assert origin == "conservative_fallback" and bs is None and npert == 1
    stability = iter([False, False, True])                 # 3 attempts: imag,imag,ok
    _s, bs, _as, _ao, origin, npert = run(b3)
    assert origin == "hessian" and calls["reopt"] == 2 and npert == 2
    assert steps == [sigma.DEFAULT_PERTURB_STEP_A,         # escalating kick
                     2 * sigma.DEFAULT_PERTURB_STEP_A], steps
    assert opre == ["qmi_reopt1", "qmi_reopt2"], opre      # distinct opt namespaces
    assert hpre == ["qmi_hess", "qmi_reopt1_hess", "qmi_reopt2_hess"], hpre
  finally:
    (qm.run_hessian, qm.run_constrained_opt, sigma.read_xtb_hessian,
     sigma.seminario_sigmas, sigma.imaginary_mode_check,
     sigma.imaginary_mode_perturbation) = saved


def exercise_provenance():
  """Provenance comment blocks (plan Sec 9.3): global (versions injected, no xTB
  subprocess) and per-site."""
  g = provenance.global_block(
    method="GFN2-xTB", solvent="ether", temperature=298.15, model_name="1abc",
    protonation="user-supplied (protonation off)",
    xtb_ver="6.7.1", cctbx_ver="2026.7")
  assert any("GFN2-xTB" in ln and "ether" in ln for ln in g)
  assert any("xtb 6.7.1" in ln and "cctbx 2026.7" in ln for ln in g)
  assert any("model: 1abc" in ln for ln in g)
  assert any(ln.startswith("protonation:") for ln in g)
  assert not any("resolution" in ln for ln in g)          # omitted when None

  st = provenance.settings_block(
    engine="orca", radius=4.0, free_radius=5.0,
    freeze_waters="peripheral", freeze_far_carboxylate_o=False, optimise_all=False,
    optimise_hydrogens="free_atoms", optimisation="looser",
    hessian="xtb", basis_set=None, nproc=6)
  assert any(ln.startswith("region:") and "radius=4.0" in ln
             and "free_radius=5.0" in ln and "freeze_waters=peripheral" in ln
             for ln in st), st
  assert not any("depth" in ln for ln in st), st
  assert any(ln.startswith("engine: orca") and "H=free_atoms" in ln
             and "opt=looser" in ln and "hessian=xtb" in ln
             and "nproc=6" in ln for ln in st), st
  assert not any("convergence=" in ln for ln in st)       # collapsed into opt
  assert not any("basis=" in ln for ln in st)             # omitted when None
  stx = provenance.settings_block(                        # xTB: no ORCA knobs
    engine="xtb", radius=4.0, free_radius=5.0,
    freeze_waters="peripheral", freeze_far_carboxylate_o=False, optimise_all=False,
    nproc=1)
  assert any(ln.startswith("engine: xtb") and "all non-cap H free" in ln
             and "nproc=1" in ln for ln in stx), stx

  bvs_est = group_args(best=3, valences={}, n_used=4)
  s = "\n".join(provenance.site_block(
    "Fe", "Fe_A54", 3, 6, 1, "guessed", "medium",
    "coordinating tyrosinate", "high-spin", bvs_est, "hessian", 1, 2))
  assert "site Fe_A54: Fe ox=+3 mult=6 charge=+1" in s
  assert "source=guessed/medium" in s
  assert "oxidation (coordinating tyrosinate); spin (high-spin)" in s
  assert "BVS most-likely ox=+3" in s
  assert "sigmas: hessian; recovered after 2 perturbation(s)" in s
  assert "1 symmetry donor(s) skipped" in s
  # conservative fallback, no perturb / no reasons / no BVS -> compact
  s2 = "\n".join(provenance.site_block(
    "Zn", "Zn_A262", 2, 1, 2, "guessed", "high", None, None, None,
    "conservative_fallback", 0, 0))
  assert "sigmas: conservative_fallback" in s2
  assert "perturbation" not in s2 and "BVS" not in s2


# ---------------------------------------------------------------------------
# qm: optimise mask + frozen/free partition
# ---------------------------------------------------------------------------
def exercise_optimise_mask():
  assert qm.optimise_mask_from_frozen(5, {1, 3}) == \
      [True, False, True, False, True]
  assert qm.optimise_mask_from_frozen(3, set()) == [True, True, True]


def exercise_sanity_bounds():
  # windows are ordered; None for an out-of-scope metal or an untabulated donor.
  tw = sanity.target_window("Zn", "N")
  assert tw is not None and tw[0] < tw[1], tw
  bw = sanity.bvs_window("Zn", "N")
  assert bw is not None and bw[0] < bw[1], bw
  assert sanity.target_window("Zn", "C") is None       # no C target
  assert sanity.bvs_window("Ca", "O") is None          # Ca not a predictor prior
  assert sanity.check_bond("Ca", "O", 2.4) == []       # out-of-scope: no crash

  # a typical bond passes; a dissociated / collapsed one trips both windows.
  assert sanity.check_bond("Zn", "N", 2.05) == []
  assert set(sanity.check_bond("Zn", "O", 3.5)) == {"atypical", "implausible"}
  assert set(sanity.check_bond("Zn", "S", 1.5)) == {"atypical", "implausible"}
  assert sanity.check_bond("Cu", "S", 2.5) == []       # real Cu-Met(S) distance

  r = sanity.check_site("Zn", [("N", 2.05), ("O", 3.5), ("S", 1.5)])
  assert r.n_atypical == 2 and r.n_implausible == 2, (r.n_atypical, r.n_implausible)
  assert len(r.details) == 2

  # provenance emits the flags when present, and nothing when clean.
  flagged = provenance.site_block(
    "Zn", "Zn_A1", 2, 1, 0, "guessed", "high", "ox", "spin", None,
    "hessian", 0, 0, sanity=r)
  assert any("sanity 2 atypical + 2 implausible" in ln for ln in flagged)
  assert any(ln.strip().startswith("sanity:") for ln in flagged)
  clean = sanity.check_site("Zn", [("N", 2.05), ("S", 2.31)])
  assert not any("sanity" in ln for ln in provenance.site_block(
    "Zn", "Zn_A1", 2, 1, 0, "guessed", "high", "ox", "spin", None,
    "hessian", 0, 0, sanity=clean))


def exercise_determine_frozen_atoms():
  model = _model(_FREEZE_PDB)
  frozen, free = qm.determine_frozen_atoms(
    model, seed_iseq=0, cap_iseqs=(), free_radius=5.0, log=null_out())
  # metal (0), near N (1), near H (3) free; far water O (2) frozen.
  assert frozen == {2}, frozen
  assert free == {0, 1, 3}, free

  # constraints.txt inventory: count header + one FROZEN/free line per atom.
  path = os.path.join(tempfile.mkdtemp(), "constraints.txt")
  qm.write_constraints_summary(model, frozen, free, path)
  lines = open(path).read().splitlines()
  assert lines[0] == "# QM optimization constraint list"
  assert "# Total atoms: 4" in lines
  assert "# Free:   3" in lines and "# Frozen: 1" in lines
  body = [ln for ln in lines if ln[:4].strip().isdigit()]
  assert len(body) == 4
  assert body[2].split()[1] == "FROZEN"                  # far water O
  assert all(body[i].split()[1] == "free" for i in (0, 1, 3))

  # cap anchors (from endo_exo) freeze the heavy atom each cap is bonded to.
  # Treat HE2 (iseq 3) as a cap anchored to NE2 (iseq 1): both end up frozen.
  frozen_c, free_c = qm.determine_frozen_atoms(
    model, seed_iseq=0, cap_iseqs=(3,), cap_anchor_iseqs=(1,),
    free_radius=5.0, log=null_out())
  assert 3 in frozen_c and 1 in frozen_c, frozen_c       # cap + its anchor
  assert free_c == {0}, free_c                           # only the metal free

  # freeze_waters modes: 'peripheral' (default) frees a coordinating water O but
  # freezes a non-coordinating one inside the relax zone; 'all' freezes both;
  # 'none' freezes neither. A water beyond free_radius is always frozen by the
  # scaffold rule regardless of mode.
  w = _model(_WATER_MODE_PDB)
  COORD, ANNULUS, BEYOND = 1, 2, 3        # ~2.1 A, ~4 A (within), ~7 A (beyond)
  for mode, coord_f, annulus_f in (("none", False, False),
                                   ("peripheral", False, True),
                                   ("all", True, True)):
    fr, _fe = qm.determine_frozen_atoms(
      w, seed_iseq=0, cap_iseqs=(), free_radius=5.0, freeze_waters=mode,
      log=null_out())
    assert (COORD in fr) == coord_f, (mode, "coord", fr)
    assert (ANNULUS in fr) == annulus_f, (mode, "annulus", fr)
    assert BEYOND in fr, (mode, "beyond frozen by radius", fr)
  # overlap fix: under 'all' the water rule freezes COORD + ANNULUS (2), NOT the
  # beyond-radius water (it is outside the free set the rule iterates).
  buf = StringIO()
  qm.determine_frozen_atoms(w, seed_iseq=0, cap_iseqs=(), free_radius=5.0,
                            freeze_waters="all", log=buf)
  assert "water O frozen (freeze_waters=all): 2" in buf.getvalue(), buf.getvalue()


# ---------------------------------------------------------------------------
# restraints: detection, assembly, .edits
# ---------------------------------------------------------------------------
def exercise_detect_first_shell():
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  donors = restraints.detect_first_shell(model, 0, sites_cart=sites)
  # 3 N (iseq 1,2,3) + water O (4); NOT the ring C (5) or the H (6).
  assert donors == [1, 2, 3, 4], donors

  # A per-element cutoff callable is honoured: admit N, reject O -> the water O
  # drops out.  (Proves detect_first_shell asks the callable per donor element.)
  only_n = lambda m_el, d_el: 5.0 if d_el.strip().upper() == 'N' else 0.0
  assert restraints.detect_first_shell(
    model, 0, sites_cart=sites, coord_cutoff=only_n) == [1, 2, 3]
  # A real Cordero cutoff tightened to alpha=1.0 excludes every contact here
  # (each N/O sits beyond r_Zn + r_donor).
  tight = util.make_coordination_cutoff('cordero', 1.0)
  assert restraints.detect_first_shell(
    model, 0, sites_cart=sites, coord_cutoff=tight) == []


def exercise_coordination_cutoff():
  fixed = util.make_coordination_cutoff()
  cordero = util.make_coordination_cutoff('cordero', 1.3)
  # 'fixed' is element-blind: the flat DEFAULT_COORDINATION_CUTOFF_A.
  assert fixed('Fe', 'O') == util.DEFAULT_COORDINATION_CUTOFF_A == 3.0
  assert fixed('Zn', 'S') == 3.0
  # 'cordero' is alpha*(r_metal + r_donor); Fe/Mn/Co take the high-spin radius,
  # the rest the cctbx.eltbx (Cordero low-spin / default) value.
  assert approx_equal(cordero('Fe', 'O'), 1.3 * (1.52 + 0.66), eps=1e-6)
  assert approx_equal(cordero('Zn', 'S'), 1.3 * (1.22 + 1.05), eps=1e-6)
  assert util.covalent_radius('Fe') == 1.52       # high-spin override
  assert util.covalent_radius('Mn') == 1.61
  assert approx_equal(util.covalent_radius(' o '), 0.66, eps=1e-2)  # eltbx, any case


def exercise_metal_ligand_restraints():
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  donors = restraints.detect_first_shell(model, 0, sites_cart=sites)

  rr = restraints.metal_ligand_restraints(model, sites, 0, donors,
                                          log=null_out())
  assert len(rr.bonds) == 4
  assert len(rr.angles) == 6           # C(4, 2)
  assert rr.n_skipped_symmetry == 0
  b0 = rr.bonds[0]
  assert approx_equal(b0["distance_ideal"], 2.0, eps=1e-3), b0   # His94 at 2.0
  # no Seminario sigmas -> conservative M-L fallback (looser than covalent 0.02/3)
  assert restraints.DEFAULT_BOND_SIGMA == 0.05
  assert restraints.DEFAULT_ANGLE_SIGMA == 5.0
  assert b0["sigma"] == restraints.DEFAULT_BOND_SIGMA
  assert b0["sigma_origin"] == "conservative_fallback"
  assert rr.angles[0]["sigma"] == restraints.DEFAULT_ANGLE_SIGMA
  assert rr.angles[0]["sigma_origin"] == "conservative_fallback"
  assert "resname ZN" in b0["atom_selection_1"]
  assert "name NE2" in b0["atom_selection_2"]

  # per-restraint (Seminario) sigmas override the conservative fallback
  bond_sigmas = dict.fromkeys(donors, 0.033)
  angle_sigmas = {frozenset((donors[i], donors[j])): 4.2
                  for i in range(len(donors))
                  for j in range(i + 1, len(donors))}
  rr2 = restraints.metal_ligand_restraints(
    model, sites, 0, donors, bond_sigmas=bond_sigmas,
    angle_sigmas=angle_sigmas, log=null_out())
  assert all(b["sigma"] == 0.033 and b["sigma_origin"] == "hessian"
             for b in rr2.bonds)
  assert all(a["sigma"] == 4.2 and a["sigma_origin"] == "hessian"
             for a in rr2.angles)
  # His94(+x)-Zn-His96(+y) is a right angle
  assert any(approx_equal(a["angle_ideal"], 90.0, eps=1e-2)
             for a in rr2.angles)

  # a contaminated angle (key absent from angle_sigmas) is still emitted, at the
  # conservative default: emit_angles never drops an angle
  partial = dict(angle_sigmas)
  partial.pop(frozenset((donors[0], donors[1])))
  rr3 = restraints.metal_ligand_restraints(
    model, sites, 0, donors, angle_sigmas=partial,
    angle_origins={k: "hessian" for k in partial}, log=null_out())
  assert len(rr3.angles) == 6          # all six emitted, none omitted
  fallback = [a for a in rr3.angles
              if a["sigma"] == restraints.DEFAULT_ANGLE_SIGMA]
  assert len(fallback) == 1 and fallback[0]["sigma_origin"] == "conservative_fallback"

  # emit_angles=False -> bonds only, no angle restraints
  rr4 = restraints.metal_ligand_restraints(
    model, sites, 0, donors, angle_sigmas=angle_sigmas, emit_angles=False,
    log=null_out())
  assert rr4.angles == [] and len(rr4.bonds) == 4

  # emitted PHIL is syntactically valid and carries the expected blocks
  edits = restraints.as_edits_phil(rr2, header="unit test")
  libtbx.phil.parse(edits)             # raises if malformed
  assert edits.count("bond {") == 4
  assert edits.count("angle {") == 6
  assert "distance_ideal" in edits and "angle_ideal" in edits


def exercise_symmetry_donor_skip():
  # A donor whose chain is absent from the original model (a symmetry image on a
  # renamed chain) has no plain input-model selection.
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  atoms = model.get_hierarchy().atoms()
  donors = restraints.detect_first_shell(model, 0, sites_cart=sites)
  keep = restraints.metal_ligand_restraints(
    model, sites, 0, donors, original_chain_ids=["A"], log=null_out())
  assert keep.n_skipped_symmetry == 0 and len(keep.bonds) == 4
  # Without provenance a symmetry-image donor is skipped.
  skip = restraints.metal_ligand_restraints(
    model, sites, 0, donors, original_chain_ids=["Q"], log=null_out())
  assert skip.n_skipped_symmetry == 4 and len(skip.bonds) == 0
  # With provenance it is emitted as a bond to its ASU parent via
  # symmetry_operation (no angle -- angles cannot cross symmetry). Provenance
  # carries the ASU parent's identity tuple, not the atom object.
  def _ident(a):
    ag = a.parent(); rg = ag.parent()
    return (rg.parent().id.strip(), rg.resseq.strip(), ag.resname.strip(),
            a.name.strip(), ag.altloc.strip())
  prov = {d: (_ident(atoms[d]), "-x,-y,-z") for d in donors}
  sym = restraints.metal_ligand_restraints(
    model, sites, 0, donors, original_chain_ids=["Q"],
    sym_image_provenance=prov, log=null_out())
  assert sym.n_skipped_symmetry == 0 and len(sym.bonds) == 4 and not sym.angles
  assert all(b["symmetry_operation"] == "-x,-y,-z" for b in sym.bonds)
  edits = restraints.as_edits_phil(sym)
  libtbx.phil.parse(edits)              # symmetry_operation renders as valid PHIL
  assert edits.count("symmetry_operation = -x,-y,-z") == 4


# ---------------------------------------------------------------------------
# protonation: metal-donor deprotonation (Cys/Tyr/Asp/Glu strip + His tautomer)
# ---------------------------------------------------------------------------
def exercise_deprotonation():
  # HIP His (both ring N protonated) coordinating via NE2 + a coordinating Cys.
  hier = _model(_HIS_HIP_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(hier, log=null_out()) == (1, 1)
  names = {a.name.strip() for a in hier.atoms()}
  assert "HE2" not in names            # coordinating-N proton removed
  assert "HD1" in names                # partner-N proton kept (neutral His)
  assert "HG" not in names             # Cys thiol stripped to anionic

  # wrong tautomer: proton only on the coordinating N -> relocated to partner.
  h2 = _model(_HIS_TAUT_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(h2, log=null_out()) == (0, 1)
  names2 = {a.name.strip() for a in h2.atoms()}
  assert "HE2" not in names2 and "HD1" in names2   # HE2 relocated + renamed

  # no metals -> no-op, proton untouched
  h3 = _model(_NOMETAL_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(h3, log=null_out()) == (0, 0)
  assert "HE2" in {a.name.strip() for a in h3.atoms()}


def exercise_deprotonation_donor_families():
  # bidentate Asp: the shared HD2 is stripped exactly once (dedup), 1 donor.
  hier = _model(_ASP_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(hier, log=null_out()) == (1, 0)
  assert "HD2" not in {a.name.strip() for a in hier.atoms()}
  # deuterated Cys thiol: DG is stripped like HG.
  h2 = _model(_CYS_D_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(h2, log=null_out()) == (1, 0)
  assert "DG" not in {a.name.strip() for a in h2.atoms()}
  # metal-scope: with metals={"Fe"} the Zn donor is out of scope -> no-op.
  h4 = _model(_ASP_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(
    h4, metals={"Fe"}, log=null_out()) == (0, 0)
  assert "HD2" in {a.name.strip() for a in h4.atoms()}


def exercise_water_protonation():
  hier = _model(_WATER_METAL_PDB).get_hierarchy()
  assert protonation.protonate_waters(hier, log=null_out()) == 2
  assert protonation.protonate_waters(hier, log=null_out()) == 0   # idempotent
  atoms = list(hier.atoms())
  o = next(a for a in atoms if a.name.strip() == "O")
  zn = next(a for a in atoms if a.element.strip().capitalize() == "Zn")
  ov = (matrix.col(zn.xyz) - matrix.col(o.xyz)).normalize()
  # cation-aware: both H sit in the hemisphere away from the metal.
  for hh in [a for a in atoms if a.element.strip().upper() == "H"]:
    hv = (matrix.col(hh.xyz) - matrix.col(o.xyz)).normalize()
    assert hv.dot(ov) <= 0.0, hv.dot(ov)


# ---------------------------------------------------------------------------
# Program: PHIL + validate (the full run needs xtb and is exercised manually)
# ---------------------------------------------------------------------------
def exercise_program():
  master = libtbx.phil.parse(master_phil_str)
  p = master.extract()
  assert p.qm_region.radius == 4.0 and p.qm.solvent == "ether"
  assert not hasattr(p.qm_region, "max_search_depth")
  assert p.qm.temperature == 100 and p.qm.nproc == 1
  assert p.qm.perturb_attempts == 1
  assert p.qm.freeze_waters == "peripheral"
  assert p.qm.angle_restraints == "auto"
  assert p.qm.work_dir == "qm_work_dir"
  assert p.output.write_files is True
  assert p.output.combine is False and p.protonation.enable is False
  assert p.protonation.water.refine.max_sweeps == 5
  assert p.protonation.water.basin.rounds == 0
  assert p.predict_only is False and p.states == 'predicted'
  assert p.use_bvs is False

  # no model -> Sorry
  prog = Program(DataManager(["model"]), p, master_phil=master,
                 logger=null_out())
  try:
    prog.validate(); raise AssertionError("expected Sorry on missing model")
  except Sorry:
    pass

  # invalid selection -> Sorry (raised before the xtb-availability check)
  dm = DataManager(["model"])
  dm.add_model("m", _model(_ZN_SITE_PDB)); dm.set_default_model("m")
  p2 = master.extract(); p2.sites.selection = ["not a valid ((( selection"]
  prog2 = Program(dm, p2, master_phil=master, logger=null_out())
  try:
    prog2.validate(); raise AssertionError("expected Sorry on bad selection")
  except Sorry:
    pass

  # unprotonated model + protonate=False -> Sorry (before the xtb check)
  dm3 = DataManager(["model"])
  dm3.add_model("noh", _model(_NOH_PDB)); dm3.set_default_model("noh")
  prog3 = Program(dm3, master.extract(), master_phil=master, logger=null_out())
  try:
    prog3.validate(); raise AssertionError("expected Sorry on unprotonated")
  except Sorry as e:
    assert "protonat" in str(e).lower(), str(e)

  # predict_only + a multi-state scan is contradictory -> Sorry.
  p4 = master.extract(); p4.predict_only = True; p4.states = 'all'
  prog4 = Program(dm, p4, master_phil=master, logger=null_out())
  try:
    prog4.validate(); raise AssertionError("expected Sorry on predict_only+scan")
  except Sorry as e:
    assert "predict_only" in str(e), str(e)

  # Bare-O water: has_hd() passes (His H) but a water has no hydrogens, so the
  # guard must still abort under protonation.enable=False. _ZN_SITE_PDB's HOH is
  # bare; adding one H (a hydroxo ligand) is accepted.
  prog5 = Program(dm, master.extract(), master_phil=master, logger=null_out())
  assert prog5._has_unprotonated_water(_model(_ZN_SITE_PDB)) is True
  try:
    prog5._require_hydrogens(_model(_ZN_SITE_PDB))
    raise AssertionError("expected Sorry on a bare-O water")
  except Sorry as e:
    assert "water" in str(e).lower(), str(e)
  oh = _ZN_SITE_PDB.replace(
    "HETATM    5  O   HOH A 263      -1.400  -1.400   0.000  1.00 10.00           O",
    "HETATM    5  O   HOH A 263      -1.400  -1.400   0.000  1.00 10.00           O\n"
    "HETATM    8  H1  HOH A 263      -2.100  -1.700   0.000  1.00 10.00           H")
  assert prog5._has_unprotonated_water(_model(oh)) is False   # hydroxo, one H
  prog5._require_hydrogens(_model(oh))                         # no raise
  # predict_only is QM-free -> validate() skips the hydrogen/water guard entirely
  # (dm holds the bare-water _ZN_SITE_PDB).
  p6 = master.extract(); p6.predict_only = True
  Program(dm, p6, master_phil=master, logger=null_out()).validate()   # no raise


def exercise_state_override():
  """User metal_state overrides: authoritative, per-selection, with fallback to
  the electronic-engine guess for unset fields / unmatched sites."""
  m = _model(_ZN_SITE_PDB)          # ZN is atom iseq 0
  def ms(**k):
    d = dict(selection=None, oxidation_state=None, multiplicity=None)
    d.update(k)
    return group_args(**d)
  b = QMRestraintsBuilder(logger=null_out())

  # No override: the predictor commits. Single-state metals -> singlet (Zn/Mg).
  st = b._resolve_state("Zn", None, "tetrahedral", 4, None)
  assert st.oxidation_state == 2 and st.multiplicity == 1 and not st.skipped
  assert st.state_source == "guessed" and st.confidence == "high", st.confidence
  assert b._resolve_state("Mg", None, "octahedral", 6, None).multiplicity == 1
  # Cu narrows by geometry: linear -> Cu(I) singlet; octahedral -> Cu(II) doublet.
  assert b._resolve_state("Cu", None, "linear", 2, None).multiplicity == 1
  cu2 = b._resolve_state("Cu", None, "octahedral", 6, None)
  assert cu2.oxidation_state == 2 and cu2.multiplicity == 2
  # Mn commits the dominant Mn(II) default at low confidence (no discriminator).
  mn = b._resolve_state("Mn", None, "octahedral", 6, None)
  assert mn.oxidation_state == 2 and mn.multiplicity == 6 and mn.confidence == "low"
  # Fe is now COMMITTED, never skipped: with no donor summary the ferrous lean
  # commits Fe(II) high-spin at low confidence.
  fe = b._resolve_state("Fe", None, "octahedral", 6, None)
  assert (fe.oxidation_state == 2 and fe.multiplicity == 5 and not fe.skipped
          and fe.confidence == "low")
  # The only skip left: a metal outside the priors with no override (Ni).
  assert b._resolve_state("Ni", None, "octahedral", 6, None).skipped

  # Full override: authoritative; runs a metal the predictor cannot guess (Ni),
  # and overrides Fe to a specific state.
  st = b._resolve_state("Ni", ms(oxidation_state=3, multiplicity=4),
                        "octahedral", 6, None)
  assert st.oxidation_state == 3 and st.multiplicity == 4 and not st.skipped
  assert st.state_source == "user-supplied" and st.confidence == "high"
  fe = b._resolve_state("Fe", ms(oxidation_state=3, multiplicity=6),
                        "octahedral", 6, None)
  assert fe.oxidation_state == 3 and fe.multiplicity == 6 and not fe.skipped

  # Partial override (oxidation only): spin is the ligand field for THAT
  # oxidation. Mn(III) high-spin octahedral d4 -> 4 unpaired -> mult 5.
  st = b._resolve_state("Mn", ms(oxidation_state=3), "octahedral", 6, None)
  assert st.oxidation_state == 3 and st.multiplicity == 5, st.multiplicity
  assert st.state_source == "user-supplied"

  # The lookup-override machinery (selection matching + conflict) is unchanged.

  # A selection that misses the metal -> no override.
  b4 = QMRestraintsBuilder(
    metal_states=[ms(selection="chain B", oxidation_state=3, multiplicity=4)],
    logger=null_out())
  assert b4._lookup_override(m, 0) is None

  # Two blocks both matching the same metal -> Sorry (disjoint required).
  b5 = QMRestraintsBuilder(
    metal_states=[ms(selection="element Zn", oxidation_state=2, multiplicity=1),
                  ms(selection="chain A and resseq 262", oxidation_state=3,
                     multiplicity=4)], logger=null_out())
  try:
    b5._lookup_override(m, 0); raise AssertionError("expected Sorry on conflict")
  except Sorry:
    pass

  # use_bvs: the BVS most-likely oxidation replaces the predictor's and the spin
  # is re-derived for it. Off (default), the predictor still wins.
  bvs3 = group_args(best=3, valences={2: 2.6, 3: 3.1}, n_used=4)
  assert b._resolve_state("Fe", None, "octahedral", 6, None,
                          bvs_est=bvs3).oxidation_state == 2   # use_bvs off
  bb = QMRestraintsBuilder(use_bvs=True, logger=null_out())
  st = bb._resolve_state("Fe", None, "octahedral", 6, None, bvs_est=bvs3)
  assert st.oxidation_state == 3 and st.multiplicity == 6      # Fe(III) HS d5
  assert st.state_source == "bvs"                             # distinct tag
  assert st.oxidation_reason == "BVS most-likely" and st.confidence == "medium"
  # A user oxidation still wins over BVS; no BVS estimate falls back to predictor.
  st = bb._resolve_state("Fe", ms(oxidation_state=2), "octahedral", 6, None,
                         bvs_est=bvs3)
  assert st.oxidation_state == 2 and st.state_source == "user-supplied"
  assert bb._resolve_state("Fe", None, "octahedral", 6, None,
                           bvs_est=None).oxidation_state == 2


def exercise_shape():
  """CShM classifier: self-consistency, a published cross-value anchor, and the
  angular-only (bond-length-invariant) property. Pure geometry, no model."""
  # Every ideal polyhedron scores ~0 against itself.
  for name, verts in shape.IDEAL_POLYHEDRA.items():
    fit = shape.csm_against_shape(verts, name)
    assert fit.csm < 1e-6, (name, fit.csm)
  # Published anchor: a perfect tetrahedron measured against square-planar is
  # exactly 100/3 (Pinsky-Avnir / SHAPE), canonical centroid-normalised.
  td = shape.IDEAL_POLYHEDRA["tetrahedral"]
  s = shape.csm_against_shape(td, "square_planar", centroid_normalized=True).csm
  assert abs(s - 100.0 / 3.0) < 1e-3, s
  # A perfect octahedron is octahedral; and against the equilateral trigonal
  # prism it reproduces the published cross-value S = 16.737 (SHAPE TPR-6).
  best = shape.find_best_shape(shape.IDEAL_POLYHEDRA["octahedral"])
  assert best.shape == "octahedral" and best.csm < 1e-6
  assert abs(best.all_fits["trigonal_prism"].csm - 16.737) < 0.05, \
    best.all_fits["trigonal_prism"].csm

  # The symmetry-reduced track (opt-in) gives the same result as full N!.
  # octahedron proper-rotation group order 24 -> 720/24 = 30 labelings.
  assert len(shape._reduced_perms("octahedral")) * 24 == 720
  noisy = flex.vec3_double(
    [(matrix.col(v) + matrix.col((0.1 * (k % 3 - 1), -0.08 * (k % 2), 0.05 * k))
      ).elems for k, v in enumerate(shape.IDEAL_POLYHEDRA["octahedral"])])
  full = shape.find_best_shape(noisy, symmetry_reduced=False)
  red = shape.find_best_shape(noisy, symmetry_reduced=True)
  assert red.shape == full.shape and abs(red.csm - full.csm) < 1e-9, \
    (full.csm, red.csm)
  # Angular-only: stretch one tetrahedral bond 3x; after unit-normalisation the
  # class is unchanged and the CShM stays ~0.
  verts = list(td)
  verts[0] = tuple(3.0 * c for c in verts[0])
  units = flex.vec3_double(
    [(matrix.col(v) / matrix.col(v).length()).elems for v in verts])
  fit = shape.find_best_shape(units)
  assert fit.shape == "tetrahedral" and fit.csm < 1e-6, (fit.shape, fit.csm)


def exercise_shape_carboxylate():
  """Carboxylate handling: a symmetric chelate collapses to one vertex, a
  monodentate one does not, an anisobidentate one drops its far O, and a
  C-terminal Asp/Glu chelating via its backbone carboxylate collapses too."""
  # Look donors up by name (mmtbx.model sorts atoms, so iseqs are not input
  # order).
  def iseqs(atoms, *names):
    want = set(names)
    return [i for i, a in enumerate(atoms) if a.name.strip() in want]

  atoms = _model(_CARBOX_SHELL_PDB).get_hierarchy().atoms()
  metal = iseqs(atoms, "ZN")[0]
  od1, od2 = iseqs(atoms, "OD1")[0], iseqs(atoms, "OD2")[0]
  ne2 = iseqs(atoms, "NE2")           # the two His N
  # Symmetric bidentate Asp (both OD) + 2 N -> raw 4, effective 3, one bidentate.
  res = shape.classify_geometry(atoms, metal, [od1, od2] + ne2)
  assert res.raw_cn == 4 and res.cn == 3, (res.raw_cn, res.cn)
  assert res.n_bidentate_carboxylate == 1
  assert "bidentate_carboxylate" in res.flags
  assert res.shape is not None
  # Monodentate Asp (only OD1) + 2 N -> no regularisation, CN 3.
  res2 = shape.classify_geometry(atoms, metal, [od1] + ne2)
  assert res2.raw_cn == 3 and res2.cn == 3
  assert res2.n_bidentate_carboxylate == 0

  # Anisobidentate Asp (2.1 / 2.7 A) + 2 N -> far O dropped, not a chelate.
  aa = _model(_CARBOX_ANISO_PDB).get_hierarchy().atoms()
  am = iseqs(aa, "ZN")[0]
  aod1, aod2 = iseqs(aa, "OD1")[0], iseqs(aa, "OD2")[0]
  ane = iseqs(aa, "NE2")
  r3 = shape.classify_geometry(aa, am, [aod1, aod2] + ane)
  assert r3.raw_cn == 4 and r3.cn == 3, (r3.raw_cn, r3.cn)
  assert r3.n_bidentate_carboxylate == 0
  assert "carboxylate_asymmetric" in r3.flags

  # C-terminal Asp chelating via its BACKBONE carboxylate (O/OXT) collapses like a
  # sidechain one, though keying the pair on resname alone would force OD1/OD2.
  ct = _model(_CTERM_CARBOX_PDB).get_hierarchy().atoms()
  cm = iseqs(ct, "ZN")[0]
  co, coxt = iseqs(ct, "O")[0], iseqs(ct, "OXT")[0]
  cne = iseqs(ct, "NE2")
  r4 = shape.classify_geometry(ct, cm, [co, coxt] + cne)
  assert r4.raw_cn == 4 and r4.cn == 3, (r4.raw_cn, r4.cn)
  assert r4.n_bidentate_carboxylate == 1
  assert "bidentate_carboxylate" in r4.flags


def exercise_electronic():
  """Ligand-field spin rules (derived, no data file): d-count, unpaired
  electrons per geometry/regime, multiplicity, and the two predicates."""
  # d-count = group - oxidation.
  assert el.d_electron_count("Zn", 2) == 10
  assert el.d_electron_count("Mn", 2) == 5
  assert el.d_electron_count("Fe", 3) == 5
  assert el.d_electron_count("Cu", 2) == 9
  assert el.d_electron_count("Mg", 2) == 0
  # Octahedral high-spin / low-spin filling across d0..d10 (the reference
  # ligand-field unpaired counts).
  oct_hs = [el.unpaired_electrons("octahedral", d) for d in range(11)]
  oct_ls = [el.unpaired_electrons("octahedral", d, low_spin=True)
            for d in range(11)]
  assert oct_hs == [0, 1, 2, 3, 4, 5, 4, 3, 2, 1, 0], oct_hs
  assert oct_ls == [0, 1, 2, 3, 2, 1, 0, 1, 2, 1, 0], oct_ls
  # Tetrahedral is always high spin (low_spin flag ignored).
  assert el.unpaired_electrons("tetrahedral", 6, low_spin=True) == 4
  # Square-planar: only strong-field d8 (diamagnetic) is a defined value.
  assert el.unpaired_electrons("square_planar", 8) == 0
  # Multiplicity = unpaired + 1: Mn(II) d5 HS sextet, Zn(II) d10 singlet,
  # Fe(II) d6 low-spin singlet.
  assert el.multiplicity(el.unpaired_electrons("octahedral", 5)) == 6
  assert el.multiplicity(el.unpaired_electrons("octahedral", 10)) == 1
  assert el.multiplicity(
    el.unpaired_electrons("octahedral", 6, low_spin=True)) == 1
  # Field-dependence: only octahedral d4-d7 is a real HS/LS choice.
  assert el.is_field_dependent("octahedral", 5)
  assert not el.is_field_dependent("octahedral", 3)
  assert not el.is_field_dependent("octahedral", 8)
  assert not el.is_field_dependent("tetrahedral", 5)
  # Encoded coverage: oct/tet all d; square-planar d8 (diamagnetic) + d9 (Cu(II)
  # doublet); field-independent d-counts (0-3, 9, 10) at any geometry; the
  # field-dependent d4-d7 at other geometries fall back (not encoded).
  assert el.is_encoded("octahedral", 7) and el.is_encoded("tetrahedral", 4)
  assert el.is_encoded("square_planar", 8) and el.is_encoded("square_planar", 9)
  assert not el.is_encoded("square_planar", 5)
  assert el.is_encoded("trigonal_bipyramidal", 10)     # field-independent d10
  assert not el.is_encoded("trigonal_bipyramidal", 5)  # field-dependent -> fallback
  # Convenience state ties it together and matches the former hardcoded states.
  st = el.spin_state("Mn", 2, "octahedral")   # high-spin default
  assert st.d_count == 5 and st.multiplicity == 6
  assert st.field_dependent and st.encoded
  assert el.spin_state("Zn", 2, "tetrahedral").multiplicity == 1

  # Oxidation-state priors (empirical, cited): accessible states + default.
  assert el.accessible_oxidation_states("Zn") == (2,)
  assert el.accessible_oxidation_states("Mn") == (2, 3)
  assert el.accessible_oxidation_states("Fe") == (2, 3)
  assert el.accessible_oxidation_states("Cu") == (1, 2)
  assert el.accessible_oxidation_states("Ni") is None   # not in scope
  assert el.default_oxidation("Zn") == 2 and el.default_oxidation("Cu") == 2
  assert el.default_oxidation("Ni") is None
  # Priors + spin rules together reproduce the former hardcoded states for the
  # in-scope metals (default oxidation, high-spin): Zn -> singlet, Mn -> sextet.
  for elem, geom, mult in (("Zn", "tetrahedral", 1), ("Mn", "octahedral", 6),
                           ("Mg", "octahedral", 1)):
    ox = el.default_oxidation(elem)
    assert el.spin_state(elem, ox, geom).multiplicity == mult, elem

  # Case/whitespace normalisation.
  assert el.default_oxidation(" zn ") == 2
  assert el.d_electron_count(" zn ", 2) == 10
  # Fe now defaults to the ferrous (photoreduction) lean, not None.
  assert el.default_oxidation("Fe") == 2
  # Error paths.
  for bad in (lambda: el.d_electron_count("Xx", 2),        # unknown element
              lambda: el.d_electron_count("Zn", -1),       # d13, out of range
              lambda: el.unpaired_electrons("octahedral", 11)):  # d11
    try:
      bad(); raise AssertionError("expected KeyError/ValueError")
    except (KeyError, ValueError):
      pass

  # ---- committing predictor (plan Sec 6): always commits, robust features only.
  def summ(spec):
    """Build a donors.summarize() result from (element, resname, name, res_id)."""
    ds = []
    for elx, rn, nm, rid in spec:
      c = D.classify(elx, rn)
      ds.append(group_args(element=c.element, kind=c.kind,
                           charge=c.charge, strong_axial=c.strong_axial,
                           res_id=rid))
    return D.summarize(ds)

  # predict_oxidation: geometry/CN commit Cu; donor hardness commits Fe.
  assert el.predict_oxidation("Zn", "tetrahedral", 4).oxidation == 2
  assert el.predict_oxidation("Cu", "linear", 2).oxidation == 1
  assert el.predict_oxidation("Cu", "trigonal", 3).oxidation == 1
  assert el.predict_oxidation("Cu", "octahedral", 6).oxidation == 2
  cu4 = el.predict_oxidation("Cu", "tetrahedral", 4)
  assert cu4.oxidation == 2 and cu4.confidence == el.CONF_LOW    # CN4 ambiguous
  assert el.predict_oxidation("Ni", "octahedral", 6) is None     # no prior
  # Fe donor-hardness commits (never bond lengths):
  tyr = summ([("O", "TYR", "OH", ("A", 1)), ("N", "HIS", "NE2", ("A", 2)),
              ("O", "ASP", "OD1", ("A", 3))])
  assert el.predict_oxidation("Fe", "octahedral", 6, tyr).oxidation == 3  # tyrosinate
  triad = summ([("N", "HIS", "NE2", ("A", 1)), ("N", "HIS", "NE2", ("A", 2)),
                ("O", "ASP", "OD1", ("A", 3))])
  ft = el.predict_oxidation("Fe", "octahedral", 6, triad)
  assert ft.oxidation == 2 and ft.confidence == el.CONF_HIGH     # facial triad
  # Rubredoxin-type Fe(S-Cys)4 is FERRIC as isolated [76, 77, 78]. This assertion
  # previously read ox == 2, from a soft-S guard whose own citation [49] never
  # supported the ferrous direction; see exercise_rubredoxin_ferric for the full
  # rule and the donor sets it must NOT catch.
  rubr = summ([("S", "CYS", "SG", ("A", i)) for i in (1, 2, 3, 4)])
  fr = el.predict_oxidation("Fe", "tetrahedral", 4, rubr)
  assert fr.oxidation == 3 and fr.confidence == el.CONF_LOW      # rubredoxin-type
  fb = el.predict_oxidation("Fe", "octahedral", 6, None)
  assert fb.oxidation == 2 and fb.confidence == el.CONF_LOW      # ferrous lean

  # Heme oxidation (Sec 6.3): the axial ligand commits ferrous/ferric. Ferric on
  # any diagnostic ferric axial ligand (phenolate/thiolate/aqua), else a LOW
  # ferrous lean (single-His/bis-His/His-Met, indeterminate).
  def heme_ax(*extra):
    return summ([("N", "HEM", "NA", ("H", 1)), ("N", "HEM", "NB", ("H", 1)),
                 ("N", "HEM", "NC", ("H", 1)), ("N", "HEM", "ND", ("H", 1))]
                + list(extra))
  h_tyr = heme_ax(("O", "TYR", "OH", ("A", 1)))                 # catalase-type
  o_tyr = el.predict_oxidation("Fe", "octahedral", 5, h_tyr)
  assert o_tyr.oxidation == 3 and o_tyr.confidence == el.CONF_MEDIUM
  h_cys = heme_ax(("S", "CYS", "SG", ("A", 1)))                 # P450-type thiolate
  o_cys = el.predict_oxidation("Fe", "octahedral", 5, h_cys)
  assert o_cys.oxidation == 3 and o_cys.confidence == el.CONF_MEDIUM
  h_aq = heme_ax(("N", "HIS", "NE2", ("A", 1)), ("O", "HOH", "O", ("W", 1)))
  o_aq = el.predict_oxidation("Fe", "octahedral", 6, h_aq)     # aquomet
  assert o_aq.oxidation == 3 and o_aq.confidence == el.CONF_MEDIUM
  h_deoxy = heme_ax(("N", "HIS", "NE2", ("A", 1)))              # 5-coordinate His
  o_dx = el.predict_oxidation("Fe", "octahedral", 5, h_deoxy)
  assert o_dx.oxidation == 2 and o_dx.confidence == el.CONF_LOW    # indeterminate lean
  h_bis = heme_ax(("N", "HIS", "NE2", ("A", 1)), ("N", "HIS", "NE2", ("A", 2)))
  assert el.predict_oxidation("Fe", "octahedral", 6, h_bis).oxidation == 2  # carrier lean

  # Mn oxidation (Sec 6.3): Mn(II) default; the SOD-fold 5-coordinate TBP
  # (3-His/carboxylate/solvent) -> Mn(III) MEDIUM; a coordinating tyrosinate ->
  # Mn(III) LOW. A 6-coordinate cupin look-alike (3-His-1-carboxylate) stays Mn(II).
  sod = summ([("N", "HIS", "NE2", ("A", 1)), ("N", "HIS", "NE2", ("A", 2)),
              ("N", "HIS", "NE2", ("A", 3)), ("O", "ASP", "OD1", ("A", 4)),
              ("O", "HOH", "O", ("W", 1))])
  m_sod = el.predict_oxidation("Mn", el.TRIGONAL_BIPYRAMIDAL, 5, sod)
  assert m_sod.oxidation == 3 and m_sod.confidence == el.CONF_MEDIUM
  cupin = summ([("N", "HIS", "NE2", ("A", 1)), ("N", "HIS", "NE2", ("A", 2)),
                ("N", "HIS", "NE2", ("A", 3)), ("O", "GLU", "OE1", ("A", 4)),
                ("O", "HOH", "O", ("W", 1)), ("O", "HOH", "O", ("W", 2))])
  assert el.predict_oxidation("Mn", "octahedral", 6, cupin).oxidation == 2  # cupin -> Mn(II)
  m_tyr = summ([("O", "TYR", "OH", ("A", 1)), ("N", "HIS", "NE2", ("A", 2)),
                ("O", "ASP", "OD1", ("A", 3))])
  m_ty = el.predict_oxidation("Mn", "octahedral", 6, m_tyr)
  assert m_ty.oxidation == 3 and m_ty.confidence == el.CONF_LOW
  m_def = el.predict_oxidation("Mn", "octahedral", 6, None)
  assert m_def.oxidation == 2 and m_def.confidence == el.CONF_LOW
  # Mn spin: always high-spin, Mn(II) d5 sextet(6), Mn(III) d4 quintet(5).
  assert el.predict_state("Mn", "octahedral", 6, cupin).multiplicity == 6
  assert el.predict_state("Mn", el.TRIGONAL_BIPYRAMIDAL, 5, sod).multiplicity == 5
  assert el.predict_spin("Mn", 3, "octahedral", None).low_spin is False

  # predict_spin: non-heme Fe high-spin; heme axial-gated (bis-His LS, aquomet HS).
  assert el.predict_spin("Fe", 2, "octahedral", triad).low_spin is False
  heme = summ([("N", "HEM", "NA", ("H", 1)), ("N", "HEM", "NB", ("H", 1)),
               ("N", "HEM", "NC", ("H", 1)), ("N", "HEM", "ND", ("H", 1)),
               ("N", "HIS", "NE2", ("A", 1)), ("N", "HIS", "NE2", ("A", 2))])
  assert heme.heme and heme.n_axial == 2 and heme.axial_all_strong
  assert el.predict_spin("Fe", 3, "octahedral", heme).low_spin is True   # bis-His
  heme_aq = summ([("N", "HEM", "NA", ("H", 1)), ("N", "HEM", "NB", ("H", 1)),
                  ("N", "HEM", "NC", ("H", 1)), ("N", "HEM", "ND", ("H", 1)),
                  ("N", "HIS", "NE2", ("A", 1)), ("O", "HOH", "O", ("W", 1))])
  assert el.predict_spin("Fe", 3, "octahedral", heme_aq).low_spin is False  # aquomet
  # His/Met is the other strong axial pair -> low spin (thioether S is strong).
  heme_hm = summ([("N", "HEM", "NA", ("H", 1)), ("N", "HEM", "NB", ("H", 1)),
                  ("N", "HEM", "NC", ("H", 1)), ("N", "HEM", "ND", ("H", 1)),
                  ("N", "HIS", "NE2", ("A", 1)), ("S", "MET", "SD", ("A", 2))])
  assert heme_hm.n_axial == 2 and heme_hm.axial_all_strong
  assert el.predict_spin("Fe", 3, "octahedral", heme_hm).low_spin is True   # His/Met
  # 5-coordinate heme (single axial) stays high spin.
  heme_5c = summ([("N", "HEM", "NA", ("H", 1)), ("N", "HEM", "NB", ("H", 1)),
                  ("N", "HEM", "NC", ("H", 1)), ("N", "HEM", "ND", ("H", 1)),
                  ("N", "HIS", "NE2", ("A", 1))])
  assert heme_5c.heme and heme_5c.n_axial == 1
  assert el.predict_spin("Fe", 3, "octahedral", heme_5c).low_spin is False  # 5-coord

  # predict_state ties it together and ALWAYS commits (None only off-scope).
  pz = el.predict_state("Zn", "tetrahedral", 4)
  assert pz.oxidation_state == 2 and pz.multiplicity == 1
  assert pz.confidence == el.CONF_HIGH and pz.low_spin is False
  pm = el.predict_state("Mn", "octahedral", 6)
  assert (pm.oxidation_state == 2 and pm.multiplicity == 6
          and pm.confidence == el.CONF_LOW)          # Mn(II) default, low conf
  # Fe is committed now (no skip): non-heme facial triad -> Fe(II) high-spin.
  pf = el.predict_state("Fe", "octahedral", 6, triad)
  assert pf.oxidation_state == 2 and pf.multiplicity == 5 and pf.low_spin is False
  # heme bis-His -> low-spin ferrous-lean d6 singlet.
  ph = el.predict_state("Fe", "octahedral", 6, heme)
  assert ph.multiplicity == 1 and ph.low_spin is True
  assert el.predict_state("Ni", "octahedral", 6) is None       # off-scope

  # enumerate_states scan enumerator. geometry=None (states=all): geometry-blind,
  # both spin extremes per accessible oxidation.
  assert el.enumerate_states("Zn") == [(2, 1)]                  # d10 singlet
  assert sorted(el.enumerate_states("Fe")) == [(2, 1), (2, 5), (3, 2), (3, 6)]
  assert el.enumerate_states("Cu") == [(1, 1), (2, 2)]         # d10 singlet, d9 doublet
  assert el.enumerate_states("Ni") is None                     # no prior
  # geometry given (states=reasonable): only the ligand-field-allowed spin(s).
  # Tetrahedral Fe is always high-spin, so the low-spin extremes are pruned.
  assert el.enumerate_states("Fe", "tetrahedral") == [(2, 5), (3, 6)]
  # Octahedral d5/d6 are genuinely field-ambiguous -> the same 4 states as blind.
  assert sorted(el.enumerate_states("Fe", "octahedral")) == \
      [(2, 1), (2, 5), (3, 2), (3, 6)]
  # An unclassified geometry cannot prune -> degrades to the blind scan.
  assert el.enumerate_states("Fe", None) == el.enumerate_states("Fe")


def exercise_deprotonation_metal_gating():
  """Cys/Tyr deprotonation is gated to anion-inducing metals, carboxylate is
  universal. Near a hard Mg the Cys stays protonated but the Asp deprotonates."""
  hier = _model(_MG_DONORS_PDB).get_hierarchy()
  protonation.deprotonate_metal_donors(hier, metals=None)   # all metals in scope
  names = [a.name.strip() for a in hier.atoms()]
  assert "HG" in names, "Cys near hard Mg must stay protonated"
  assert "HD2" not in names, "Asp carboxylate must deprotonate (metal-agnostic)"


def exercise_donors():
  """Biological donor classification + coordination summary (robust features
  only: residue/element identity, no bond lengths)."""
  # classify: kind / formal charge / strong-axial role.
  assert D.classify("N", "HIS").kind == "his_n"
  assert D.classify("N", "HIS").strong_axial                 # imidazole axial
  assert D.classify("O", "ASP").charge == -1                 # carboxylate anion
  assert D.classify("O", "TYR").kind == "tyrosinate_o"
  cys = D.classify("S", "CYS")
  assert cys.charge == -1 and not cys.strong_axial           # thiolate, soft
  assert D.classify("S", "MET").strong_axial                 # thioether His/Met
  assert D.classify("O", "HOH").kind == "water_o"
  assert D.classify("N", "HEM").kind == "het_n"              # porphyrin N (het)

  def desc(elx, rn, nm, rid):
    c = D.classify(elx, rn)
    return group_args(element=c.element, kind=c.kind,
                      charge=c.charge, strong_axial=c.strong_axial, res_id=rid)

  # Heme by FeN4 topology: 4 N from ONE residue -> heme; 4 His N from four
  # residues -> NOT heme (the discriminator is same-residue, not element count).
  heme = D.summarize([desc("N", "HEM", nm, ("H", 1))
                      for nm in ("NA", "NB", "NC", "ND")]
                     + [desc("N", "HIS", "NE2", ("A", 1)),
                        desc("N", "HIS", "NE2", ("A", 2))])
  assert heme.heme and heme.n_axial == 2 and heme.axial_all_strong
  four_his = D.summarize([desc("N", "HIS", "NE2", ("A", i)) for i in range(4)])
  assert not four_his.heme
  # Motif counts + hard-anion charge dedup per residue (bidentate Asp = one -1).
  s = D.summarize([desc("O", "ASP", "OD1", ("A", 5)),
                   desc("O", "ASP", "OD2", ("A", 5)),   # same Asp -> one group
                   desc("O", "TYR", "OH", ("A", 6)),
                   desc("N", "HIS", "NE2", ("A", 7))])
  assert s.n_carboxylate == 1 and s.n_tyrosinate == 1 and s.n_his == 1
  assert s.q_hard == -2                                  # carboxylate + tyrosinate


def exercise_predict_only():
  """predict_only: detection + committing prediction, stopping before the QM
  step (no xtb, no restraints emitted)."""
  res = QMRestraintsBuilder(predict_only=True, logger=null_out()).run(
    _model(_ZN_SITE_PDB))
  assert len(res) == 1
  r = res[0]
  assert not r.skipped and r.restraints is None and r.edits is None
  assert r.oxidation_state == 2 and r.multiplicity == 1
  assert r.confidence == "high" and r.charge is None
  assert r.site_id == "ZN_chainA_res262", r.site_id


def exercise_bvs():
  """Bond-valence-sum oxidation estimate: CIF parse, valence sum, and the
  deduce-most-self-consistent-state estimator with its missing-parameter guard."""
  # Parse: Zn(II)-O picks the Brown-Altermatt row; Fe(II)/Fe(III) R0 differ.
  assert bvs.valence_params("Zn", 2, "O") == (1.704, 0.37)
  assert bvs.valence_params("Fe", 2, "O")[0] != bvs.valence_params("Fe", 3, "O")[0]
  # A tetrahedral Zn(II)-O4 at ~1.97 A sums to V ~ 2.
  V, n = bvs.valence_sum("Zn", 2, [("O", 1.97)] * 4)
  assert n == 4 and abs(V - 2.0) < 0.1, (V, n)
  # estimate_oxidation: Zn only has +2 accessible; V ~ 2 so best is +2.
  e = bvs.estimate_oxidation("Zn", (2,), [("O", 1.97)] * 4)
  assert e is not None and e.best == 2 and e.n_used == 4
  assert abs(e.valences[2] - 2.0) < 0.1, e.valences
  # Fe-O6 short bonds (~2.00 A octahedral) are more self-consistent with +3:
  # V_III lands nearer 3 than V_II lands to 2, so the estimator deduces +3.
  fe = bvs.estimate_oxidation("Fe", (2, 3), [("O", 2.00)] * 6)
  assert fe is not None and set(fe.valences) == {2, 3} and fe.n_used == 6
  assert fe.best == 3, fe.valences
  # Longer Fe-O6 bonds (~2.15 A) tip the estimate the other way, to +2.
  assert bvs.estimate_oxidation("Fe", (2, 3), [("O", 2.15)] * 6).best == 2
  # No vendored parameters for any candidate -> None; valence_sum guard.
  assert bvs.estimate_oxidation("Zn", (2,), [("Xx", 2.0)]) is None
  assert bvs.valence_sum("Zn", 2, [("Xx", 2.0)]) == (None, 0)


def exercise_bvs_hygiene():
  """BVS reliability guards: unchecked-row selection, the close-call margin, and
  the asymmetric-coverage refusal (with the None-verdict paths that refusal
  opens up)."""
  # bvparm2020.cif lists rows in decreasing reliability, but its header warns
  # that rows marked 'unchecked' are unreliable. Where a checked alternative
  # follows, it wins. Both of these are load-bearing for qmi.
  assert bvs.valence_params("Fe", 2, "S") == (2.125, 0.37)   # not the 2.12 unchecked row
  assert bvs.valence_params("Cu", 1, "O") == (1.601, 0.335)  # not (1.610, 0.37); B changes too

  # Margin = the gap in |V - q| between the winner and the runner-up. A
  # well-separated site is not a close call; a marginal one is. The threshold
  # matters because BVS discriminates on ~0.01-0.06 A of bond length, less than
  # the spread between re-refinements of one crystal [84].
  wide = bvs.estimate_oxidation("Fe", (2, 3), [("O", 2.15)] * 6)
  assert wide.best == 2 and wide.margin > bvs.CLOSE_CALL_MARGIN
  assert not wide.close_call, wide.margin
  tight = bvs.estimate_oxidation("Fe", (2, 3), [("O", 2.09)] * 6)
  assert tight.margin < bvs.CLOSE_CALL_MARGIN and tight.close_call, tight.margin

  # A single accessible state has nothing to compare against, so there is no
  # margin and the result must not masquerade as a discriminated verdict.
  one = bvs.estimate_oxidation("Zn", (2,), [("O", 1.97)] * 4)
  assert one.best == 2 and one.margin is None and not one.close_call

  # Asymmetric parameter coverage: Mn(2+) has an S row, Mn(3+) does not, so a
  # mixed Mn thiolate site would sum the two candidates over different bond
  # counts and silently push the verdict to Mn(II). Refuse instead.
  mn = bvs.estimate_oxidation("Mn", (2, 3), [("S", 2.40), ("N", 2.20), ("O", 2.10)])
  assert mn.coverage_asymmetric and mn.best is None
  assert mn.n_by_state[2] != mn.n_by_state[3], mn.n_by_state

  # best=None must not crash its consumers. Provenance reports the refusal, and
  # _resolve_state falls through to the predictor instead of committing None,
  # with the source label following the branch that actually decided.
  lines = provenance.site_block("Mn", "Mn_A1", 2, 6, 0, "guessed", "low",
                                "r", "s", mn, "hessian", 0, 0)
  assert any("no verdict" in ln for ln in lines), lines
  b = QMRestraintsBuilder(use_bvs=True, logger=null_out())
  st = b._resolve_state("Mn", None, "octahedral", 6, None, bvs_est=mn)
  assert st.oxidation_state is not None and st.state_source == "guessed"

  # The margin and the close-call warning reach the header.
  note = provenance.site_block("Fe", "Fe_A1", 2, 5, 0, "guessed", "low", "r",
                               "s", tight, "hessian", 0, 0)
  assert any("close call" in ln for ln in note), note


def exercise_unclassified_donors():
  """Donors the residue-name classifier could not place must be counted and
  surfaced, and must lower confidence ONLY where the oxidation actually came
  from the donor summary."""
  def desc(elx, rn, rid):
    c = D.classify(elx, rn)
    return group_args(element=c.element, kind=c.kind, charge=c.charge,
                      strong_axial=c.strong_axial, res_id=("A", rid, "", rn))

  # A ligand no whitelist covers (citrate here, but the rule is general: any
  # het_* kind) contributes to no motif count, so the summary carries nothing.
  cit = D.summarize([desc("O", "FLC", "402")] * 3 + [desc("O", "FLC", "403")] * 3)
  assert cit.n_unclassified == 6 and not cit.motif_evidence
  assert cit.unclassified == ["FLC.O"], cit.unclassified

  # A heme's four pyrrole N ARE explained, by the FeN4 topology test, so they
  # are not reported as unclassified. But `heme` is NOT motif evidence: it is
  # derived from unplaced het_n donors, so a 100%-het macrocycle must still read
  # as uninformative.
  heme = D.summarize([desc("N", "HEM", "350")] * 4
                     + [desc("N", "HIS", "169"), desc("O", "HOH", "678")])
  assert heme.heme and heme.n_unclassified == 0 and heme.motif_evidence
  chl = D.summarize([desc("N", "CLA", "1")] * 4)
  assert chl.heme and not chl.motif_evidence, "heme must not count as evidence"
  # Only the four N the test consumed are excused; a fifth N on the same residue
  # (an intramolecular axial, as in cobalamin) is still unexplained.
  cbl = D.summarize([desc("N", "B12", "1")] * 5)
  assert cbl.n_unclassified == 1, cbl.n_unclassified

  # Partial: motif evidence survives, so the site is reported but not capped.
  part = D.summarize([desc("O", "ASP", "60"), desc("O", "TYR", "92"),
                      desc("N", "HIS", "250"), desc("O", "CO3", "402")])
  assert part.n_unclassified == 1 and part.motif_evidence

  b = QMRestraintsBuilder(logger=null_out())
  # Fires: Fe reads the donor summary, and here it is empty.
  fe = b._resolve_state("Fe", None, "octahedral", 6, cit)
  assert fe.confidence == el.CONF_LOW
  assert "NO donor was classified" in fe.oxidation_reason, fe.oxidation_reason
  # Does NOT fire where the donor evidence is irrelevant to the oxidation: a
  # single accessible state (Mg/Zn) is definitional, and Cu decides on CN and
  # geometry. Capping those would contradict their own reason strings.
  mg = b._resolve_state("Mg", None, "octahedral", 6,
                        D.summarize([desc("O", "GNP", "201")] * 6))
  assert mg.confidence == el.CONF_HIGH, mg.oxidation_reason
  assert "NO donor was classified" not in (mg.oxidation_reason or "")
  cu = b._resolve_state("Cu", None, "octahedral", 6,
                        D.summarize([desc("O", "LIG", "1")] * 6))
  assert cu.confidence == el.CONF_HIGH, cu.oxidation_reason
  # The helper the cap keys on must track predict_oxidation's branches.
  assert el.oxidation_uses_donors("Fe") and el.oxidation_uses_donors("Mn")
  assert not el.oxidation_uses_donors("Zn") and not el.oxidation_uses_donors("Cu")

  # And it reaches the .edits header, so an empty summary is never silent.
  lines = provenance.site_block("Fe", "Fe_A401", 2, 5, 0, "guessed", "low",
                                fe.oxidation_reason, "s", None, "hessian", 0, 0,
                                donor_summary=cit)
  assert any("unclassified donors: 6" in ln for ln in lines), lines


def exercise_rubredoxin_ferric():
  """Rubredoxin-type Fe(S-Cys)4 commits FERRIC, and the rule is narrow enough to
  leave every neighbouring donor set alone."""
  def desc(elx, rn, rid):
    c = D.classify(elx, rn)
    return group_args(element=c.element, kind=c.kind, charge=c.charge,
                      strong_axial=c.strong_axial, res_id=("A", rid, "", rn))
  def ox(cn, ds, geom=None):
    return el.predict_oxidation(
      "Fe", geom or ("tetrahedral" if cn == 4 else "octahedral"), cn,
      D.summarize(ds))

  # Exactly four Cys thiolates, nothing else: ferric as isolated [76, 77, 78],
  # LOW confidence because the site is ferrous in reservoir roles [87] and
  # partly photoreduced at routine dose [84].
  rd = ox(4, [desc("S", "CYS", str(i)) for i in (5, 8, 38, 41)])
  assert rd.oxidation == 3, rd.reason
  assert rd.confidence == el.CONF_LOW                  # never medium or high
  assert "rubredoxin-type" in rd.reason
  assert "ferrous" in rd.reason and "photoreduced" in rd.reason  # names the ambiguity
  # High-spin d5 -> sextet, from the same donor set.
  st = el.predict_state("Fe", "tetrahedral", 4,
                        D.summarize([desc("S", "CYS", str(i)) for i in (5, 8, 38, 41)]))
  assert st.oxidation_state == 3 and st.multiplicity == 6

  # NOT caught: nitrile hydratase (3 Cys-derived S + backbone N, CN 5) is a
  # LOW-SPIN Fe(III) [86]; committing it ferric here would be right for the wrong
  # reason and would still carry the wrong spin, so the rule must not fire.
  nhase = ox(5, [desc("S", "CYS", str(i)) for i in (1, 2, 3)]
             + [desc("N", "ALA", "4"), desc("N", "ALA", "5")])
  assert "rubredoxin-type" not in nhase.reason, nhase.reason

  # NOT caught: every [2Fe-2S] iron presents 2 Cys S + 2 inorganic S, which the
  # previous >=2-Cys form of this rule swallowed. (The mononuclear gate upstream
  # is the real defence; this is belt and braces.)
  fes = ox(4, [desc("S", "CYS", "1"), desc("S", "CYS", "2"),
               desc("S", "SF4", "9"), desc("S", "SF4", "9")])
  assert "rubredoxin-type" not in fes.reason, fes.reason

  # NOT caught: a fifth donor means it is not a rubredoxin-type site.
  plus_w = ox(5, [desc("S", "CYS", str(i)) for i in (1, 2, 3, 4)]
              + [desc("O", "HOH", "9")])
  assert "rubredoxin-type" not in plus_w.reason, plus_w.reason
  # NOT caught: superoxide-reductase centre II (4 His + 1 Cys) is ferrous.
  sor = ox(5, [desc("N", "HIS", str(i)) for i in (1, 2, 3, 4)]
           + [desc("S", "CYS", "5")])
  assert sor.oxidation == 2, sor.reason

  # Unchanged neighbours: the two HIGH-confidence Fe rules must not move.
  triad = ox(6, [desc("N", "HIS", "1"), desc("N", "HIS", "2"),
                 desc("O", "ASP", "3")]
             + [desc("O", "HOH", str(i)) for i in (7, 8, 9)])
  assert triad.oxidation == 2 and triad.confidence == el.CONF_HIGH
  tyr = ox(6, [desc("O", "TYR", "92"), desc("O", "TYR", "191"),
               desc("O", "ASP", "60"), desc("N", "HIS", "250"),
               desc("O", "CO3", "402"), desc("O", "CO3", "402")])
  assert tyr.oxidation == 3 and tyr.confidence == el.CONF_HIGH


def exercise_dose_pair_invariance():
  """A radiation-damaged twin must get the SAME committed state, and the
  advisory disagreement that follows must be reported rather than buried.

  9TA4 and 9TA6 are one rubredoxin crystal measured at 50 kGy and then at 1 MGy
  [84]. Every refinement-robust feature is identical between them (4 Cys S,
  CN 4, tetrahedral, charge -1); only the Fe-S distances differ, by 0.037 A,
  which matches a full chemical reduction of the same protein (5NW3 -> 5OME,
  0.0368 A). So the predictor MUST return the same state for both: that is the
  invariant, not an oversight. Letting the distances decide instead would be
  circular, since qmi's own restraint moves them [84].
  """
  def desc(elx, rn, rid):
    c = D.classify(elx, rn)
    return group_args(element=c.element, kind=c.kind, charge=c.charge,
                      strong_axial=c.strong_axial, res_id=("A", rid, "", rn))
  cys4 = [desc("S", "CYS", str(i)) for i in (6, 9, 39, 42)]
  summary = D.summarize(cys4)

  # Same donors, same answer, deliberately: dose is not a donor property.
  low, high = (el.predict_oxidation("Fe", "tetrahedral", 4, summary)
               for _ in range(2))
  assert low.oxidation == high.oxidation == 3, (low.oxidation, high.oxidation)
  assert low.confidence == high.confidence == el.CONF_LOW
  assert low.reason == high.reason
  assert "photoreduced" in low.reason, low.reason

  # The observed geometries DO differ, and BVS sees it. Real deposited Fe-S
  # distances, measured from 9TA4 and 9TA6.
  d_low = (2.2821, 2.2573, 2.3204, 2.2620)      # 9TA4, 50 kGy,  mean 2.2804
  d_high = (2.3243, 2.2910, 2.3528, 2.2999)     # 9TA6, 1 MGy,   mean 2.3170
  est_low = bvs.estimate_oxidation("Fe", (2, 3), [("S", d) for d in d_low])
  est_high = bvs.estimate_oxidation("Fe", (2, 3), [("S", d) for d in d_high])
  assert est_low.best == 3 and not est_low.close_call, est_low.margin
  assert est_high.best == 2 and est_high.close_call, est_high.margin
  # ...and the high-dose verdict is a knife edge: an order of magnitude less
  # separation than the low-dose one, well under CLOSE_CALL_MARGIN.
  assert est_high.margin < 0.1 < est_low.margin, (est_low.margin,
                                                  est_high.margin)

  # The disagreement must reach the .edits header, with the close-call caveat
  # and the override, and must NOT appear where BVS agrees.
  hi = provenance.site_block("Fe", "Fe_A101", 3, 6, -1, "guessed", "low",
                             low.reason, "s", est_high, "hessian", 0, 0)
  assert any("BVS DISAGREES" in ln for ln in hi), hi
  assert any("close call" in ln for ln in hi), hi
  assert any("metal_state" in ln for ln in hi), hi
  lo = provenance.site_block("Fe", "Fe_A101", 3, 6, -1, "guessed", "low",
                             low.reason, "s", est_low, "hessian", 0, 0)
  assert not any("BVS DISAGREES" in ln for ln in lo), lo
  # A user who already chose the state is told of the disagreement but not
  # advised to use the override they have plainly already found.
  usr = provenance.site_block("Fe", "Fe_A101", 3, 6, -1, "user-supplied", "",
                              low.reason, "s", est_high, "hessian", 0, 0)
  assert any("BVS DISAGREES" in ln for ln in usr), usr
  assert not any("metal_state" in ln for ln in usr), usr
  # And a BVS-sourced state cannot disagree with itself.
  src = provenance.site_block("Fe", "Fe_A101", 2, 5, -2, "bvs", "low",
                              low.reason, "s", est_high, "hessian", 0, 0)
  assert not any("BVS DISAGREES" in ln for ln in src), src


def exercise_donor_perception():
  """The bond-graph perception layer: a delocalised charge must be read off the
  resonance GROUP, since a Lewis walk parks it on one arbitrary member."""
  from mmtbx.geometry_restraints.qmi import donor_perception as P
  # Acetate: C(methyl)-C(carboxyl)(-O)(-O), charge parked on ONE oxygen.
  ace = P.DonorPerception(["C", "C", "O", "O"],
                          [(0, 1), (1, 2), (1, 3)], charges={2: -1})
  assert ace.resonance_group(2) == (2, 3) == ace.resonance_group(3)
  # The point of the whole exercise: the uncharged oxygen is the coordinating
  # one about half the time, and reading the ATOM there gives 0.
  assert ace.charges.get(3, 0) == 0
  assert ace.group_charge(2) == -1 and ace.group_charge(3) == -1
  # Deduped by group: a bidentate carboxylate is ONE anionic group, not two.
  assert ace.group_charges([2, 3]) == {(2, 3): -1}

  # A metal is never a substituent: a donor bonded to one must not look
  # saturated, or every sigma-count rule downstream misfires.
  with_m = P.DonorPerception(["C", "C", "O", "O", "Fe"],
                             [(0, 1), (1, 2), (1, 3), (2, 4)], charges={2: -1},
                             metal_iseqs=[4])
  assert with_m.heavy_neighbours(2) == [1]
  assert with_m.resonance_group(2) == (2, 3)
  assert with_m.n_sigma(2) == 1

  # Sigma count and H: a hydroxyl O has 2 (one C, one H), an alkoxide 1.
  roh = P.DonorPerception(["C", "O", "H"], [(0, 1), (1, 2)])
  assert roh.n_sigma(1) == 2 and roh.bonded_h(1) == [2]
  assert roh.resonance_group(1) == (1,), "a hydroxyl is its own group"

  # A group on a ring the Lewis walk never solved is reported untrusted rather
  # than handed over as a result.
  ring = P.DonorPerception(["C", "O"], [(0, 1)], charges={1: -1},
                           unresolved_iseqs=[0, 1])
  assert not ring.group_is_resolved(1)
  assert ace.group_is_resolved(2), "a resolved group stays usable"
  # The PARENT counts as much as the group. A terminal O hanging off an unsolved
  # aromatic ring is its own group and never intersects the ring, so checking
  # the group alone would call an unsolved phenolate trustworthy.
  phenol = P.DonorPerception(["C", "O"], [(0, 1)], charges={1: -1},
                             unresolved_iseqs=[0])
  assert phenol.resonance_group(1) == (1,)
  assert not phenol.group_is_resolved(1), "parent on an unsolved ring"

  # The charge of a group is not always all on its terminal atoms. A nitro
  # group is formally neutral: N is +1 and one O is -1, so reading the oxygens
  # alone gives -1 for a fragment that carries no charge.
  nitro = P.DonorPerception(["C", "N", "O", "O"],
                            [(0, 1), (1, 2), (1, 3)], charges={1: 1, 2: -1})
  assert nitro.group_charge(2) == 0, nitro.group_charge(2)
  # Including the parent leaves the ordinary anions exactly where they were,
  # since their parent atoms are formally neutral...
  assert ace.group_charge(2) == -1
  pho = P.DonorPerception(["P", "O", "O", "O", "O"],
                          [(0, 1), (0, 2), (0, 3), (0, 4)],
                          charges={1: -1, 2: -1, 3: -1})
  assert pho.group_charge(1) == -3, pho.group_charge(1)
  # ...and corrects a sulfate written zwitterionically, S(2+) with four O(-),
  # which the terminal atoms alone would read as -4.
  sul = P.DonorPerception(["S", "O", "O", "O", "O"],
                          [(0, 1), (0, 2), (0, 3), (0, 4)],
                          charges={0: 2, 1: -1, 2: -1, 3: -1, 4: -1})
  assert sul.group_charge(1) == -2, sul.group_charge(1)


def exercise_bond_graph_donor_charge():
  """Anionic charge the residue names cannot see is recovered from the bond
  graph, and ONLY where the names ran out."""
  def desc(elx, rn, rid, gid=None, gq=None):
    c = D.classify(elx, rn)
    return group_args(element=c.element, kind=c.kind, charge=c.charge,
                      strong_axial=c.strong_axial, res_id=("A", rid, "", rn),
                      group_id=gid, group_charge=gq)

  # Bis-citrate Fe (6G7Q): four carboxylate groups plus two neutral hydroxyls,
  # none of them on any residue-name whitelist.
  cit = D.summarize(
    [desc("O", "FLC", "402", gid=(1,), gq=-1),
     desc("O", "FLC", "402", gid=(2,), gq=-1),
     desc("O", "FLC", "402", gid=(3,), gq=0),      # hydroxyl, neutral
     desc("O", "FLC", "403", gid=(4,), gq=-1),
     desc("O", "FLC", "403", gid=(5,), gq=-1),
     desc("O", "FLC", "403", gid=(6,), gq=0)])
  assert cit.n_anionic_unplaced == 4, cit.n_anionic_unplaced
  assert cit.q_hard == -4, cit.q_hard
  # ...and that is enough to commit the ferric state the whitelist alone missed.
  fe = el.predict_oxidation("Fe", "octahedral", 6, cit)
  assert fe.oxidation == 3 and fe.confidence == el.CONF_MEDIUM, fe.reason
  assert "bond graph" in fe.reason, fe.reason

  # A group is counted ONCE however many of its atoms coordinate: a bidentate
  # carboxylate presents -1, not -2.
  bi = D.summarize([desc("O", "FLC", "402", gid=(1,), gq=-1),
                    desc("O", "FLC", "402", gid=(1,), gq=-1)])
  assert bi.n_anionic_unplaced == 1 and bi.q_hard == -1

  # Capped at -1 per group, so the calibrated q_hard thresholds keep their
  # meaning: a carbonate at -2 does not count double.
  carb = D.summarize([desc("O", "CO3", "402", gid=(1,), gq=-2)])
  assert carb.q_hard == -1, carb.q_hard

  # Donors the residue names DID place are untouched, so nothing the whitelist
  # already gets right can be double-counted or overridden.
  asp = D.summarize([desc("O", "ASP", "60", gid=(1,), gq=-1),
                     desc("O", "TYR", "92", gid=(2,), gq=-1)])
  assert asp.n_anionic_unplaced == 0 and asp.q_hard == -2

  # Hard O only. Porphyrin pyrrole N carry the macrocycle charge, and letting
  # them into q_hard would make EVERY heme ferric; a thiolate S is soft and has
  # no business in a rule cited to hard-oxygen HSAB.
  # A 5-coordinate His heme (deoxy globin / resting peroxidase) is the sharp
  # case: it has no anionic axial ligand at all, so if the macrocycle's own
  # charge leaked into q_hard the ferric limb would fire on the ONE thing the
  # heme rules are written to ignore.
  heme = D.summarize([desc("N", "HEM", "350", gid=(i,), gq=-1)
                      for i in range(4)]
                     + [desc("N", "HIS", "169")])
  assert heme.heme and heme.n_anionic_unplaced == 0 and heme.q_hard == 0
  fe_heme = el.predict_oxidation("Fe", "square_pyramidal", 5, heme)
  assert fe_heme.oxidation == 2, fe_heme.reason
  soft = D.summarize([desc("S", "SF4", "1", gid=(i,), gq=-1) for i in range(4)])
  assert soft.n_anionic_unplaced == 0 and soft.q_hard == 0

  # No bond graph (an uninterpretable ligand): the predictor falls back to the
  # residue names it used before, and says the state is a default.
  none = D.summarize([desc("O", "FLC", "402") for _ in range(6)])
  assert none.n_anionic_unplaced == 0 and none.q_hard == 0
  assert not none.motif_evidence

  # A group id that is falsy but present must not fall back to object identity,
  # or two donors of ONE group would count as two.
  same = D.summarize([desc("O", "FLC", "402", gid=(), gq=-1),
                      desc("O", "FLC", "402", gid=(), gq=-1)])
  assert same.n_anionic_unplaced == 1, same.n_anionic_unplaced

  # The heme limb fires on a SINGLE anionic group, so a bond-graph charge is
  # load-bearing there on its own and the reason must disclose it rather than
  # asserting a phenolate it has not seen.
  heme_ox = D.summarize([desc("N", "HEM", "350", gid=(i,)) for i in range(4)]
                        + [desc("O", "LIG", "9", gid=(9,), gq=-1)])
  assert heme_ox.heme and heme_ox.q_hard == -1
  fe_ax = el.predict_oxidation("Fe", "octahedral", 6, heme_ox)
  assert fe_ax.oxidation == 3 and "bond graph" in fe_ax.reason, fe_ax.reason

  # Mn reads the donor summary but never q_hard, so a recovered charge is not
  # evidence for a Mn site and must not lift its confidence cap.
  b = QMRestraintsBuilder(logger=null_out())
  mn = b._resolve_state("Mn", None, "octahedral", 6, cit)
  assert mn.confidence == el.CONF_LOW, mn.oxidation_reason
  assert "NO donor was classified" in mn.oxidation_reason, mn.oxidation_reason
  # ...whereas for Fe, which does read it, the cap correctly lifts.
  fe_state = b._resolve_state("Fe", None, "octahedral", 6, cit)
  assert "NO donor was classified" not in (fe_state.oxidation_reason or "")
  assert el.oxidation_uses_donor_charge("Fe")
  assert not el.oxidation_uses_donor_charge("Mn")

  # The recovered charge reaches the .edits header, and replaces the
  # "state is a default" note, which would now be untrue.
  lines = provenance.site_block("Fe", "Fe_A401", 3, 6, 0, "guessed", "medium",
                                fe.reason, "s", None, "hessian", 0, 0,
                                donor_summary=cit)
  assert any("recovered from the bond graph" in ln for ln in lines), lines
  assert not any("state is a default" in ln for ln in lines), lines


def exercise_atcun_deprotonation():
  # ATCUN Cu N-terminus: the coordinating backbone amide (ALA N-H) is stripped to
  # an amidate, AND the N-terminal amine (GLY -NH3+) loses the H nearest the
  # metal (-> -NH2). Both count as donors (returns 2, 0).
  h = _model(_ATCUN_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(h, log=null_out()) == (2, 0)
  ala = {a.name.strip() for a in h.atoms()
         if a.parent().resname.strip() == "ALA"}
  gly = {a.name.strip() for a in h.atoms()
         if a.parent().resname.strip() == "GLY"}
  assert "H" not in ala, ala                             # amide -> amidate
  assert len(gly & {"H1", "H2", "H3"}) == 2, gly         # -NH3+ -> -NH2
  # metal gating: the same backbone N near a Zn is neither amidate nor amine.
  hz = _model(_ZN_AMIDE_PDB).get_hierarchy()
  assert protonation.deprotonate_metal_donors(hz, log=null_out()) == (0, 0)
  assert "H" in {a.name.strip() for a in hz.atoms()}


def exercise_freeze_far_carboxylate():
  # A monodentate Asp (OD1 in, OD2 beyond the cutoff): the far OD2 is frozen only
  # with the flag; without it OD2 stays free.
  mono = _model(_CARBOX_MONO_PDB)
  name = {a.i_seq: a.name.strip() for a in mono.get_hierarchy().atoms()}
  fr, _ = qm.determine_frozen_atoms(
    mono, seed_iseq=0, cap_iseqs=(), free_radius=5.0,
    freeze_far_carboxylate_o=True, log=null_out())
  fr0, _ = qm.determine_frozen_atoms(
    mono, seed_iseq=0, cap_iseqs=(), free_radius=5.0,
    freeze_far_carboxylate_o=False, log=null_out())
  assert {name[i] for i in fr} == {"OD2"}, {name[i] for i in fr}
  assert {name[i] for i in fr0} == set(), fr0
  # a genuine bidentate Asp (both OD within cutoff) has no far O -> untouched.
  bi = _model(_CARBOX_SHELL_PDB)
  bname = {a.i_seq: a.name.strip() for a in bi.get_hierarchy().atoms()}
  frb, _ = qm.determine_frozen_atoms(
    bi, seed_iseq=0, cap_iseqs=(), free_radius=5.0,
    freeze_far_carboxylate_o=True, log=null_out())
  assert not any(bname[i] in ("OD1", "OD2") for i in frb), frb


def exercise_coordination_mode_changes():
  # A detected donor (NE2) the opt pushes OUT -> 'left'; a non-donor water O the
  # opt pulls IN -> 'entered'.
  model = _model(_MODE_PDB)
  atoms = model.get_hierarchy().atoms()
  idx = {a.name.strip(): a.i_seq for a in atoms}
  opt = flex.vec3_double([a.xyz for a in atoms])
  opt[idx["NE2"]] = (4.0, 0.0, 0.0)             # donor leaves (> 3.0)
  opt[idx["O"]] = (2.5, 0.0, 0.0)               # non-donor enters (< 3.0)
  details = sanity.coordination_mode_changes(
    model, opt, metal_iseq=idx["ZN"], donor_iseqs={idx["NE2"]},
    coord_cutoff=None)
  assert any("left" in d and "NE2" in d for d in details), details
  assert any("entered" in d and "HOH" in d for d in details), details
  # unchanged geometry -> no mode change.
  same = sanity.coordination_mode_changes(
    model, flex.vec3_double([a.xyz for a in atoms]), metal_iseq=idx["ZN"],
    donor_iseqs={idx["NE2"]}, coord_cutoff=None)
  assert same == [], same


# ---------------------------------------------------------------------------
# ORCA backend wiring (no live QM: dispatch, input construction, .hess parse)
# ---------------------------------------------------------------------------

# ZN + a coordinating water (O + 2 H), a His donor N with a ring H, and a frozen
# capping carbon.  iseqs are file order: 0 ZN, 1 O, 2/3 water H, 4 N, 5 His H,
# 6 cap C.
_ORCA_HFREEZE_PDB = """\
CRYST1   60.000   60.000   60.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
HETATM    2  O   HOH A 263       2.100   0.000   0.000  1.00 10.00           O
HETATM    3  H1  HOH A 263       2.700   0.700   0.000  1.00 10.00           H
HETATM    4  H2  HOH A 263       2.700  -0.700   0.000  1.00 10.00           H
ATOM      5  NE2 HIS A  94       0.000   2.100   0.000  1.00 10.00           N
ATOM      6  HE2 HIS A  94       0.500   2.800   0.000  1.00 10.00           H
ATOM      7  CB  LEU A  50       0.000  -3.000   0.000  1.00 10.00           C
ATOM      8  HB2 LEU A  50       0.500  -3.500   0.500  1.00 10.00           H
ATOM      9  HC  ACE A 999       3.000   3.000   3.000  1.00 10.00           H
END
"""


def exercise_engine_options():
  x = qm.engine_options()
  assert x.engine == "xtb" and x.method == qm.DEFAULT_XTB_METHOD
  assert x.solvent_model == qm.DEFAULT_XTB_SOLVENT
  assert x.hessian == "xtb" and x.optimise_hydrogens == "free_atoms"
  o = qm.engine_options(engine="orca")
  assert o.engine == "orca" and o.method == qm.DEFAULT_ORCA_METHOD
  assert o.solvent_model == qm.DEFAULT_ORCA_SOLVENT
  assert o.optimisation == "LooseOpt" and o.convergence is None
  # explicit overrides win; a None basis normalises to "".
  o2 = qm.engine_options(engine="orca", method="TPSS", basis_set=None,
                         solvent_model="CPCM(water)", optimisation="Opt",
                         hessian="partial", optimise_hydrogens="all", exe="/x")
  assert o2.method == "TPSS" and o2.basis_set == "" and o2.hessian == "partial"
  assert o2.optimisation == "Opt" and o2.optimise_hydrogens == "all"
  assert o2.exe == "/x"


def exercise_read_orca_hessian():
  # ORCA writes the (3N x 3N) Hessian in 5-wide column blocks, each led by an
  # integer column-index header row.  Use a symmetric M[i][j]=i+j so no
  # symmetrisation is triggered.  2 atoms -> 6x6.
  n = 6
  M = [[i + j for j in range(n)] for i in range(n)]
  tmp = tempfile.mkdtemp()
  path = os.path.join(tmp, "probe.hess")
  with open(path, "w") as fh:
    fh.write("$hessian\n%d\n" % n)
    for c0 in range(0, n, 5):
      cols = list(range(c0, min(c0 + 5, n)))
      fh.write("    " + "".join("%18d" % c for c in cols) + "\n")
      for r in range(n):
        fh.write("%5d" % r
                 + "".join("%18.10E" % M[r][c] for c in cols) + "\n")
    fh.write("$vibrational_frequencies\n")
  h = sigma.read_orca_hessian(path, 2)
  assert h.focus() == (6, 6), h.focus()
  assert approx_equal(h[(0, 0)], 0.0)
  assert approx_equal(h[(0, 5)], 5.0)
  assert approx_equal(h[(5, 5)], 10.0)
  assert approx_equal(h[(2, 3)], 5.0) and approx_equal(h[(3, 2)], 5.0)
  try:
    sigma.read_orca_hessian(path, 3)      # wrong atom count must raise
  except ValueError:
    pass
  else:
    raise AssertionError("expected ValueError on size mismatch")


def exercise_free_atoms_hydrogens():
  # 'free_atoms' rule: a hydrogen is frozen iff its heavy atom is frozen, with
  # two exemptions (capping H stay frozen, coordinating-water H stay free).
  b = QMRestraintsBuilder(engine="orca", orca_exe="x", logger=null_out())
  model = _model(_ORCA_HFREEZE_PDB)
  ZN, WO, WH1, WH2, NE2, HE2, CB, HB2, CAP = range(9)
  # base partition: the water O (freeze_waters) and the scaffold CB are frozen
  # heavy atoms, the cap is frozen, every other non-cap H is free.
  frozen = set(b._freeze_fixed_atom_hydrogens(
    model, frozen_iseqs=[WO, CB, CAP], donor_iseqs=[WO, NE2], cap_iseqs=[CAP]))
  assert HB2 in frozen, "a hydrogen on a frozen heavy atom must be frozen"
  assert HE2 not in frozen, "a hydrogen on a free heavy atom must stay free"
  assert CAP in frozen, "a capping hydrogen must stay frozen"
  assert WH1 not in frozen and WH2 not in frozen, \
    "coordinating-water hydrogens stay free even though their O is frozen"
  assert {WO, CB} <= frozen, "frozen heavy atoms stay frozen"
  assert not ({ZN, NE2} & frozen), "free heavy atoms stay free"


def exercise_run_hessian_dispatch():
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  # The runners are now the managers' run_hessian methods; monkeypatch those.
  saved = (xtb_manager.xtb_manager.run_hessian,
           orca_manager.orca_manager.run_hessian,
           sigma.read_xtb_hessian, sigma.read_orca_hessian)
  calls = []
  try:
    def fake_xtb_run(self, *a, **k): calls.append("xtb_run"); return "xp"
    def fake_orca_run(self, *a, **k): calls.append("orca_run"); return "op"
    def fake_xtb_read(p, n): calls.append(("xtb_read", p)); return "XH"
    def fake_orca_read(p, n): calls.append(("orca_read", p)); return "OH"
    xtb_manager.xtb_manager.run_hessian = fake_xtb_run
    orca_manager.orca_manager.run_hessian = fake_orca_run
    sigma.read_xtb_hessian = fake_xtb_read
    sigma.read_orca_hessian = fake_orca_read
    # xtb engine -> xtb runner + xtb reader
    r = qm.run_hessian(model, sites, 2, 1, engine_opts=qm.engine_options())
    assert calls == ["xtb_run", ("xtb_read", "xp")] and r == "XH", calls
    # orca engine but hessian=xtb -> still the xtb Hessian
    calls[:] = []
    r = qm.run_hessian(model, sites, 2, 1,
                       engine_opts=qm.engine_options(engine="orca",
                                                     hessian="xtb"))
    assert calls == ["xtb_run", ("xtb_read", "xp")] and r == "XH", calls
    # orca engine, hessian=partial -> orca runner + orca reader
    calls[:] = []
    r = qm.run_hessian(model, sites, 2, 1, hessian_atoms=[0, 1],
                       engine_opts=qm.engine_options(engine="orca",
                                                     hessian="partial"))
    assert calls == ["orca_run", ("orca_read", "op")] and r == "OH", calls
  finally:
    (xtb_manager.xtb_manager.run_hessian,
     orca_manager.orca_manager.run_hessian,
     sigma.read_xtb_hessian, sigma.read_orca_hessian) = saved


def exercise_run_constrained_opt_dispatch():
  model = _model(_ZN_SITE_PDB)
  work = tempfile.mkdtemp()
  saved = (orca_manager.orca_manager, xtb_manager.xtb_manager)
  built = []
  class _Fake(object):
    def __init__(self, *a, **k): built.append((type(self).tag, k.get("nproc")))
    def set_ligand_atoms(self, m): self.mask = m
    def get_opt(self, **k):
      built.append((type(self).tag + "_opt", k.get("optimise_h"))); return None, "C"
  class FakeOrca(_Fake): tag = "orca"
  class FakeXtb(_Fake): tag = "xtb"
  try:
    orca_manager.orca_manager = FakeOrca
    xtb_manager.xtb_manager = FakeXtb
    _c, m1 = qm.run_constrained_opt(model, [2, 3, 4], 2, 1, nproc=4,
                                    work_dir=work, engine_opts=qm.engine_options())
    assert isinstance(m1, FakeXtb) and m1.mask[0] is True and m1.mask[2] is False
    _c2, m2 = qm.run_constrained_opt(
      model, [2, 3, 4], 2, 1, nproc=4, work_dir=work,
      engine_opts=qm.engine_options(engine="orca", optimisation="Opt", exe="/x"))
    assert isinstance(m2, FakeOrca) and m2.optimisation == "Opt"
    assert ("xtb", 4) in built and ("orca", 4) in built
    assert ("xtb_opt", False) in built and ("orca_opt", False) in built
  finally:
    (orca_manager.orca_manager, xtb_manager.xtb_manager) = saved


# ---------------------------------------------------------------------------
# metal_protonation (the rules, decided separately from applying them)
#
# protonation.py is exercised above through deprotonate_metal_donors, which
# walks a hierarchy. These cover the rule module underneath it, whose point is
# that decide() answers from facts alone, with no structure in sight, so a
# protonation program can call it from inside a placement loop before the
# hydrogen in question exists.
# ---------------------------------------------------------------------------
_MP_ZN_CYS_PDB = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  CB  CYS A   1       9.000  10.000  10.000  1.00 20.00           C
ATOM      2  SG  CYS A   1      10.400  10.000  10.000  1.00 20.00           S
ATOM      3  HG  CYS A   1      10.800  10.900  10.000  1.00 20.00           H
HETATM    4 ZN    ZN A 100      12.500  10.000  10.000  1.00 20.00          ZN
END
"""

_MP_NA_CYS_PDB = _MP_ZN_CYS_PDB.replace(
  "HETATM    4 ZN    ZN A 100", "HETATM    4 NA    NA A 100").replace(
  "          ZN", "          NA")

# A HIP His coordinating through ND1 (the site above uses NE2), so the
# coordinating ring N is the one the nearest-metal pick has to find.
_MP_HIS_HIP_PDB = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
ATOM      1  CG  HIS A   1       8.000  10.000  10.000  1.00 20.00           C
ATOM      2  ND1 HIS A   1       9.100   9.200  10.000  1.00 20.00           N
ATOM      3  CD2 HIS A   1       8.400  11.300  10.000  1.00 20.00           C
ATOM      4  CE1 HIS A   1      10.100  10.050  10.000  1.00 20.00           C
ATOM      5  NE2 HIS A   1       9.700  11.300  10.000  1.00 20.00           N
ATOM      6  HD1 HIS A   1       9.150   8.200  10.000  1.00 20.00           H
ATOM      7  HE2 HIS A   1      10.200  12.100  10.000  1.00 20.00           H
HETATM    8 ZN    ZN A 100      11.300   8.600  10.000  1.00 20.00          ZN
END
"""


def _mp_hier(pdb_str):
  return iotbx.pdb.input(source_info=None,
                         lines=pdb_str.splitlines()).construct_hierarchy()


def _mp_cutoff(_metal_element, _donor_element):
  return 3.0


def _mp_metals(hier):
  return [a for a in hier.atoms()
          if a.element.strip().capitalize() in ("Zn", "Na", "Cu", "Ni", "Fe")]


def exercise_rules_decide_is_pure():
  """The oracle answers from facts alone, with no structure in sight."""
  thiolate = mp.decide(mp.donor_facts("CYS", "SG", "S", "Zn", 2.3,
                                      h_names=["HG"]))
  assert thiolate.action == mp.REMOVE_H, thiolate.reason
  assert thiolate.h_name == "HG" and thiolate.mark_acceptor

  # Already anionic: nothing to do, and no pretending otherwise.
  done = mp.decide(mp.donor_facts("CYS", "SG", "S", "Zn", 2.3, h_names=[]))
  assert done.action == mp.KEEP, done.reason

  # HSAB gate: a hard, low-charge ion does not induce a thiolate [37, 38].
  soft_gate = mp.decide(mp.donor_facts("CYS", "SG", "S", "Na", 2.8,
                                       h_names=["HG"]))
  assert soft_gate.action == mp.KEEP and soft_gate.rule == "hsab-gate"
  # ...but a carboxylate is anionic by pH whatever the metal, so it is NOT
  # gated: it fires on the same hard ion the thiolate rule refuses.
  acid = mp.decide(mp.donor_facts("GLU", "OE1", "O", "Na", 2.4,
                                  h_names=["HE2"]))
  assert acid.action == mp.REMOVE_H, acid.reason


def exercise_rules_histidine_stays_neutral():
  """A coordinating His must never become an imidazolate."""
  # Both ring N protonated: drop the coordinating one, leaving a neutral ring.
  hip = mp.decide(mp.donor_facts("HIS", "ND1", "N", "Zn", 2.1,
                                 h_names=["HD1"], partner_h_names=["HE2"]))
  assert hip.action == mp.REMOVE_H and hip.rule == "his-hip"

  # Only the coordinating N protonated: the wrong tautomer. The proton MOVES
  # across the ring; removing it would strip the ring's last proton.
  taut = mp.decide(mp.donor_facts("HIS", "ND1", "N", "Zn", 2.1,
                                  h_names=["HD1"], partner_h_names=[]))
  assert taut.action == mp.MOVE_TO_PARTNER, taut.reason
  assert taut.new_name == "HE2", taut.new_name
  # Deuterium keeps its isotope through the rename.
  d_taut = mp.decide(mp.donor_facts("HIS", "NE2", "N", "Zn", 2.1,
                                    h_names=["DE2"], partner_h_names=[]))
  assert d_taut.new_name == "DD1", d_taut.new_name

  # Coordinating N already bare: correct neutral tautomer, leave it alone.
  ok = mp.decide(mp.donor_facts("HIS", "NE2", "N", "Zn", 2.1, h_names=[],
                                partner_h_names=["HD1"]))
  assert ok.action == mp.KEEP


def exercise_rules_atcun():
  """ATCUN N-terminal deprotonation, and its gate."""
  amine = mp.decide(mp.donor_facts("GLY", "N", "N", "Cu", 2.0,
                                   h_names=["H1", "H2", "H3"],
                                   is_standard_aa=True))
  assert amine.action == mp.REMOVE_H and amine.rule == "atcun-amine"
  # Which of the three equivalent H goes is geometry, not chemistry, so the
  # rule names the criterion instead of guessing an atom.
  assert amine.h_select == "nearest_to_metal" and amine.h_name is None

  amidate = mp.decide(mp.donor_facts("ALA", "N", "N", "Cu", 2.0,
                                     h_names=["H"], is_standard_aa=True))
  assert amidate.action == mp.REMOVE_H and amidate.rule == "atcun-amidate"

  # Only Cu/Ni do this biologically; a Zn backbone N keeps its amide H.
  zn = mp.decide(mp.donor_facts("ALA", "N", "N", "Zn", 2.0, h_names=["H"],
                                is_standard_aa=True))
  assert zn.action == mp.KEEP and zn.rule == "atcun-gate"


def exercise_rules_observed_hydrogens():
  """An OBSERVED hydrogen is data. The default refuses to delete it."""
  facts = mp.donor_facts("CYS", "SG", "S", "Zn", 2.3, h_names=["HG"],
                         h_observed=True)
  assert mp.decide(facts).action == mp.KEEP
  assert mp.decide(facts).rule == "observed-h"
  # A caller that knows the H is placed rather than measured can say so.
  loose = mp.default_policy(deprotonate_observed_h=True)
  assert mp.decide(facts, policy=loose).action == mp.REMOVE_H


def exercise_rules_policy_is_data():
  """The metal sets are choices a caller may override, not laws."""
  facts = mp.donor_facts("CYS", "SG", "S", "Pb", 2.6, h_names=["HG"])
  assert mp.decide(facts).action == mp.KEEP     # Pb is not in the default set
  extended = mp.default_policy(
    anion_inducing_metals=frozenset(mp.ANION_INDUCING_METALS | {"Pb"}))
  assert mp.decide(facts, policy=extended).action == mp.REMOVE_H


def exercise_rules_propose_mutates_nothing():
  """propose() is read-only; apply_changes() is what changes the structure."""
  hier = _mp_hier(_MP_ZN_CYS_PDB)
  before = hier.atoms().size()
  changes = mp.propose(hier, _mp_metals(hier), _mp_cutoff)
  assert len(changes) == 1, [e.rule for e in changes]
  assert hier.atoms().size() == before, "propose must not mutate"
  assert changes[0].action == mp.REMOVE_H and changes[0].rule == "anionic-donor"
  # The change identifies its hydrogen in a way that survives a re-sort.
  assert changes[0].h_id[4] == "CYS" and changes[0].h_id[5] == "HG"
  counts = mp.apply_changes(hier, changes)
  assert counts.n_removed == 1
  assert hier.atoms().size() == before - 1
  assert "HG" not in [a.name.strip() for a in hier.atoms()]

  # The same donor at a hard ion is left alone, end to end.
  na = _mp_hier(_MP_NA_CYS_PDB)
  assert mp.propose(na, _mp_metals(na), _mp_cutoff) == []


def exercise_rules_propose_his_tautomer():
  """A doubly-protonated coordinating His loses exactly one proton."""
  hier = _mp_hier(_MP_HIS_HIP_PDB)
  changes = mp.propose(hier, _mp_metals(hier), _mp_cutoff)
  assert len(changes) == 1, [(e.rule, e.h_id) for e in changes]
  assert changes[0].rule == "his-hip"
  # The metal is nearest ND1, so that is the coordinating ring N.
  assert changes[0].donor_atom.name.strip() == "ND1", changes[0].donor_atom.name
  assert changes[0].h_id[5] == "HD1"
  mp.apply_changes(hier, changes)
  names = [a.name.strip() for a in hier.atoms()]
  assert "HD1" not in names and "HE2" in names, names


def exercise_rules_apply_does_not_sort():
  """A caller may hold i_seq-indexed structures, so the applier must not
  silently reorder their atoms. qmi's own wrapper re-sorts afterwards; that is
  its choice to make, not the rule module's."""
  hier = _mp_hier(_MP_ZN_CYS_PDB)
  order_before = [a.name.strip() for a in hier.atoms()]
  changes = mp.propose(hier, _mp_metals(hier), _mp_cutoff)
  mp.apply_changes(hier, changes)
  order_after = [a.name.strip() for a in hier.atoms()]
  assert order_after == [n for n in order_before if n != "HG"], order_after


def exercise_rules_no_metals_is_a_noop():
  hier = _mp_hier(_MP_ZN_CYS_PDB)
  assert mp.propose(hier, [], _mp_cutoff) == []


def run():
  exercise_metals()
  exercise_hessian_parse()
  exercise_engine_options()
  exercise_read_orca_hessian()
  exercise_free_atoms_hydrogens()
  exercise_run_hessian_dispatch()
  exercise_run_constrained_opt_dispatch()
  exercise_seminario_bond()
  exercise_seminario_angle_guards()
  exercise_sigma_conversion()
  exercise_classify_angle_sigmas()
  exercise_emit_angles()
  exercise_coordination_number()
  exercise_imaginary_mode()
  exercise_perturb_reopt()
  exercise_provenance()
  exercise_sanity_bounds()
  exercise_coordination_mode_changes()
  exercise_optimise_mask()
  exercise_determine_frozen_atoms()
  exercise_freeze_far_carboxylate()
  exercise_detect_first_shell()
  exercise_coordination_cutoff()
  exercise_metal_ligand_restraints()
  exercise_symmetry_donor_skip()
  exercise_deprotonation()
  exercise_deprotonation_donor_families()
  exercise_deprotonation_metal_gating()
  exercise_atcun_deprotonation()
  exercise_rules_decide_is_pure()
  exercise_rules_histidine_stays_neutral()
  exercise_rules_atcun()
  exercise_rules_observed_hydrogens()
  exercise_rules_policy_is_data()
  exercise_rules_propose_mutates_nothing()
  exercise_rules_propose_his_tautomer()
  exercise_rules_apply_does_not_sort()
  exercise_rules_no_metals_is_a_noop()
  exercise_water_protonation()
  exercise_state_override()
  exercise_shape()
  exercise_shape_carboxylate()
  exercise_electronic()
  exercise_donors()
  exercise_predict_only()
  exercise_bvs()
  exercise_bvs_hygiene()
  exercise_unclassified_donors()
  exercise_rubredoxin_ferric()
  exercise_donor_perception()
  exercise_bond_graph_donor_charge()
  exercise_dose_pair_invariance()
  exercise_program()


if __name__ == "__main__":
  run()
  print("OK")
