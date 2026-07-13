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
from libtbx.test_utils import approx_equal
from libtbx.utils import null_out
from scitbx.array_family import flex

from mmtbx.geometry_restraints.qmi import qm, restraints, sigma, metals


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
  from mmtbx.geometry_restraints.qmi import protonation
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
  from mmtbx.geometry_restraints.qmi import protonation
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


# ---------------------------------------------------------------------------
# Program: PHIL + validate (the full run needs xtb and is exercised manually)
# ---------------------------------------------------------------------------
def exercise_program():
  from iotbx.data_manager import DataManager
  from libtbx.utils import Sorry
  from mmtbx.programs.qmi import Program, master_phil_str
  master = libtbx.phil.parse(master_phil_str)
  p = master.extract()
  assert p.buffer_radius == 3.0 and p.solvent == "ether"
  assert p.temperature == 298.15 and p.write_files is True
  assert p.combine is False and p.protonate is False

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
  p2 = master.extract(); p2.selection = ["not a valid ((( selection"]
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
  exercise_program()


if __name__ == "__main__":
  run()
  print("OK")
