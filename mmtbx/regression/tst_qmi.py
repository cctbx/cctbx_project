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
from mmtbx.geometry_restraints.qmi import (
  qm, restraints, sigma, metals, protonation, shape, bvs)
from mmtbx.geometry_restraints.qmi import electronic as el
from mmtbx.geometry_restraints.qmi import donors as D
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

# Zn + a bare coordinating water, for the cation-aware water-H placement.
_WATER_METAL_PDB = """\
CRYST1   30.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1 ZN    ZN A 262       0.000   0.000   0.000  1.00 10.00          ZN
HETATM    2  O   HOH A 300       2.100   0.000   0.000  1.00 10.00           O
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


# ---------------------------------------------------------------------------
# qm: optimise mask + frozen/free partition
# ---------------------------------------------------------------------------
def exercise_optimise_mask():
  assert qm.optimise_mask_from_frozen(5, {1, 3}) == \
      [True, False, True, False, True]
  assert qm.optimise_mask_from_frozen(3, set()) == [True, True, True]


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


# ---------------------------------------------------------------------------
# restraints: detection, assembly, .edits
# ---------------------------------------------------------------------------
def exercise_detect_first_shell():
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  donors = restraints.detect_first_shell(model, 0, sites_cart=sites)
  # 3 N (iseq 1,2,3) + water O (4); NOT the ring C (5) or the H (6).
  assert donors == [1, 2, 3, 4], donors


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
  assert b0["sigma"] == restraints.DEFAULT_BOND_SIGMA
  assert b0["sigma_origin"] == "placeholder"
  assert "resname ZN" in b0["atom_selection_1"]
  assert "name NE2" in b0["atom_selection_2"]

  # per-restraint (Seminario) sigmas override the placeholder
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

  # emitted PHIL is syntactically valid and carries the expected blocks
  edits = restraints.as_edits_phil(rr2, header="unit test")
  libtbx.phil.parse(edits)             # raises if malformed
  assert edits.count("bond {") == 4
  assert edits.count("angle {") == 6
  assert "distance_ideal" in edits and "angle_ideal" in edits


def exercise_symmetry_donor_skip():
  # A donor whose chain is absent from the original model (a symmetry image on
  # a renamed chain) has no input-model selection and must be skipped.
  model = _model(_ZN_SITE_PDB)
  sites = model.get_hierarchy().atoms().extract_xyz()
  donors = restraints.detect_first_shell(model, 0, sites_cart=sites)
  keep = restraints.metal_ligand_restraints(
    model, sites, 0, donors, original_chain_ids=["A"], log=null_out())
  assert keep.n_skipped_symmetry == 0 and len(keep.bonds) == 4
  skip = restraints.metal_ligand_restraints(
    model, sites, 0, donors, original_chain_ids=["Q"], log=null_out())
  assert skip.n_skipped_symmetry == 4 and len(skip.bonds) == 0


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
  assert p.qm.temperature == 298.15 and p.qm.nproc == 1
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
  monodentate one does not, and an anisobidentate one drops its far O."""
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
      c = D.classify(elx, rn, nm)
      ds.append(group_args(element=c.element, kind=c.kind, field=c.field,
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
  rubr = summ([("S", "CYS", "SG", ("A", i)) for i in (1, 2, 3, 4)])
  fr = el.predict_oxidation("Fe", "tetrahedral", 4, rubr)
  assert fr.oxidation == 2 and fr.confidence == el.CONF_LOW      # soft-S guard
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
  # classify: kind / field / formal charge / strong-axial role.
  assert D.classify("N", "HIS", "NE2").kind == "his_n"
  assert D.classify("N", "HIS", "NE2").strong_axial          # imidazole axial
  assert D.classify("O", "ASP", "OD1").charge == -1          # carboxylate anion
  assert D.classify("O", "TYR", "OH").kind == "tyrosinate_o"
  cys = D.classify("S", "CYS", "SG")
  assert cys.charge == -1 and not cys.strong_axial           # thiolate, soft
  assert D.classify("S", "MET", "SD").strong_axial           # thioether His/Met
  assert D.classify("O", "HOH", "O").kind == "water_o"
  assert D.classify("N", "HEM", "NA").kind == "het_n"        # porphyrin N (het)

  def desc(el, rn, nm, rid):
    c = D.classify(el, rn, nm)
    return group_args(element=c.element, kind=c.kind, field=c.field,
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


def run():
  exercise_metals()
  exercise_hessian_parse()
  exercise_seminario_bond()
  exercise_seminario_angle_guards()
  exercise_sigma_conversion()
  exercise_optimise_mask()
  exercise_determine_frozen_atoms()
  exercise_detect_first_shell()
  exercise_metal_ligand_restraints()
  exercise_symmetry_donor_skip()
  exercise_deprotonation()
  exercise_deprotonation_donor_families()
  exercise_deprotonation_metal_gating()
  exercise_water_protonation()
  exercise_state_override()
  exercise_shape()
  exercise_shape_carboxylate()
  exercise_electronic()
  exercise_donors()
  exercise_predict_only()
  exercise_bvs()
  exercise_program()


if __name__ == "__main__":
  run()
  print("OK")
