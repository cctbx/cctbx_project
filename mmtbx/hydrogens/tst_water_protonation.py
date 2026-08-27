"""Regression tests for ``water_protonation.place_water_hydrogens``: the
H-bond-aware water-hydrogen placer.

Minimal, self-contained cases exercise the parts that have been easy
to get wrong:

* **Acceptor-directed placement.** One HOH flanked by two acceptor O on the
  H-O-H cone.  Each placed proton must point at an acceptor (H1 -> nearer,
  H2 -> the other), clash-free and at the right O-H / H-O-H geometry.  This
  is also the regression for the ``memory_id`` atom-identity bug: with the
  water's own atoms not excluded from its clash search, the protons stopped
  pointing at the acceptors.

* **Cation repulsion.** A water coordinating an Mg whose only acceptor is
  on the metal side.  Both protons must still end up in the hemisphere
  *away* from the metal (pointing H+ at M(n+) is electrostatically wrong);
  acceptor-direction alone would pull one toward it.

* **Protonated N is a donor.** An N that already carries an H is excluded
  as an acceptor, so a water points at a nearby bare N instead.

* **Lone-pair-directed placement.** With ``lone_pair_directed`` the O-H aims
  at a carbonyl O's sp2 lone-pair lobe (~120 deg off the C=O axis) rather
  than its nucleus; off by default it aims at the nucleus.

* **Joint H1/H2 placement.** With two acceptors too acute (60 deg) for both
  H-bonds from an on-acceptor H1, ``joint`` balances the pair so the worse
  H-bond is much better than greedy's leftover.

* **Idempotency.** A water that already carries two H is left untouched.

* **Single-H water (possible hydroxide).** A water carrying one H is left
  untouched by default and reported; ``existing_h="complete"`` builds the
  partner on the *deposited* proton's cone (the regression for building it
  on a recomputed one, which gave an arbitrary H-O-H angle);
  ``existing_h="reorient"`` strips and re-places both.

* **Metal annotation.** The single-H report flags a coordinating cation
  using a first-shell cutoff, so a metal inside the looser proton-repulsion
  radius but beyond a bond is not flagged.

* **Reorient existing.** With ``existing_h="reorient"`` a protonated water
  whose H point the wrong way is stripped and re-placed toward its
  acceptor; the default leaves it alone.

* **HETATM record type.** Placed H inherit the parent O's HETATM flag.

* **Heavy cations repel.** The cation set is not limited to the first row:
  a Pt keeps both protons out of its hemisphere.

* **Refinement.** A tight cluster of bare waters clashes after the greedy
  pass; the relaxation sweeps must reduce the count of close H-H contacts.

* **A completed water stays completed.** In a clashing cluster where half
  the waters carry a deposited proton, the refinement sweeps and basin
  kicks must leave every deposited proton where it is and keep each new
  partner on its cone, rather than re-deriving the O-H1 axis.

* **Element override.** ``element="D"`` forces deuterium (named ``D1``/``D2``)
  onto an HOH water; the default places H.

* **Experiment detection.** ``_detect_neutron`` reads the experiment type
  (EXPDTA / ``_exptl.method``) and falls back to the presence of D atoms.

* **O-H length auto-selection.** With no ``oh_length`` the placer takes the
  neutron distance for a model carrying D and the X-ray distance otherwise;
  an explicit value overrides both.
"""

from __future__ import absolute_import, division, print_function

import math

import iotbx.pdb
from scitbx import matrix
from libtbx.utils import format_cpu_times

from mmtbx.hydrogens import water_protonation as wp


# One HOH (bare O) with two acceptor O placed on the 104.5 deg H-O-H cone:
# A1 at 2.6 A along +x, A2 at 2.8 A at 104.5 deg from +x (+z azimuth).  The
# placer should orient H1 -> A1 (nearer) and H2 -> A2.
_TWO_ACCEPTOR_PDB = """\
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  O   ACA D   1       7.600   5.000   5.000  1.00 10.00           O
HETATM    3  O   ACB D   2       4.299   5.000   7.711  1.00 10.00           O
END
"""

# Water O coordinating an Mg (2.5 A along +x). The ONLY acceptor sits on
# the Mg side (30 deg off the O->Mg axis), so acceptor-direction alone
# would pull an H toward the metal -- only the cation repulsion keeps both
# protons in the hemisphere away from it.
_MG_WATER_PDB = """\
HETATM    1 MG    MG A   1       2.500   0.000   0.000  1.00 10.00          MG
HETATM    2  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    3  O   ACA D   1       2.338   1.350   0.000  1.00 10.00           O
END
"""

# A water flanked by two N: a *nearer* protonated N (a donor: it carries an
# H) and a farther bare N (an acceptor). H1 must skip the donor and point
# at the bare N.
_DONOR_N_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  N   ACC D   1      -2.800   0.000   0.000  1.00 10.00           N
HETATM    3  N   DNR E   1       2.600   0.000   0.000  1.00 10.00           N
HETATM    4  H   DNR E   1       3.600   0.000   0.000  1.00 10.00           H
END
"""

# A water in the plane of an sp2 carbonyl (C=O with the C bonded to two
# more C, defining the plane), off to one side. With lone-pair-directed
# placement the O-H aims at an in-plane lobe (~120 deg from C=O); without
# it, at the O nucleus (a poorer C=O...H angle).
_CARBONYL_PDB = """\
HETATM    1  O   HOH W   1       2.000   1.500   0.000  1.00 10.00           O
HETATM    2  O   ACO A   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    3  C   ACO A   1      -1.220   0.000   0.000  1.00 10.00           C
HETATM    4  C   ACO A   1      -1.950   1.260   0.000  1.00 10.00           C
HETATM    5  C   ACO A   1      -1.950  -1.260   0.000  1.00 10.00           C
END
"""

# A water flanked by two acceptor O only 60 deg apart -- too acute for the
# 104.5 deg H-O-H to hit both from an on-acceptor H1. Greedy nails one and
# leaves the other poorly aligned; joint placement balances the two.
_TWO_ACCEPTOR_60_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   ACA D   1       2.800   0.000   0.000  1.00 10.00           O
HETATM    3  O   ACB E   1       1.400   2.424   0.000  1.00 10.00           O
END
"""

# A water that is already protonated (O + 2 H) -- must be left untouched.
_PROTONATED_WATER_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       0.000   0.957   0.000  1.00 10.00           H
HETATM    3  H2  HOH W   1       0.926  -0.239   0.000  1.00 10.00           H
END
"""

# A protonated water whose H point AWAY from the lone acceptor (+x): with
# --reorient-existing the H are stripped and re-placed toward the acceptor.
_BAD_PROTONATED_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  H1  HOH W   1      -0.984   0.000   0.000  1.00 10.00           H
HETATM    3  H2  HOH W   1       0.000  -0.984   0.000  1.00 10.00           H
HETATM    4  O   ACA D   1       2.700   0.000   0.000  1.00 10.00           O
END
"""

# A water carrying a single H (a common way of writing hydroxide), with one
# acceptor along +x. The deposited proton points along +y, well off the
# acceptor, so completing it must build the partner on *that* proton's cone.
_SINGLE_H_WATER_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       0.000   0.957   0.000  1.00 10.00           H
HETATM    3  O   ACA D   1       2.700   0.000   0.000  1.00 10.00           O
END
"""

# Single-H water coordinating an Mg at 2.1 A: the report must flag it. The
# FAR copy moves the Mg to 2.9 A, inside the proton-repulsion radius but
# beyond a first-shell bond, so it must not be flagged.
_MG_SINGLE_H_PDB = """\
HETATM    1 MG    MG A   1      -2.100   0.000   0.000  1.00 10.00          MG
HETATM    2  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    3  H1  HOH W   1       0.000   0.957   0.000  1.00 10.00           H
END
"""

_MG_SINGLE_H_FAR_PDB = _MG_SINGLE_H_PDB.replace("-2.100", "-2.900")

# Six bare waters on a tight 2.8 A grid with no other acceptors: the greedy
# pass leaves at least one close H-H contact that refinement relaxes.
_WATER_CLUSTER_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   HOH W   2       0.000   2.800   0.000  1.00 10.00           O
HETATM    3  O   HOH W   3       2.800   0.000   0.000  1.00 10.00           O
HETATM    4  O   HOH W   4       2.800   2.800   0.000  1.00 10.00           O
HETATM    5  O   HOH W   5       5.600   0.000   0.000  1.00 10.00           O
HETATM    6  O   HOH W   6       5.600   2.800   0.000  1.00 10.00           O
END
"""


# Half the waters of a clashing 2.6 A cube carry a deposited H. Completing
# them exercises refinement and basin-hopping with fixed-axis records in
# play: the deposited protons must not move and the partners must stay on
# their cones through every sweep and kick.
_MIXED_CLUSTER_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   HOH W   2       0.000   0.000   2.600  1.00 10.00           O
HETATM    3  H1  HOH W   2       0.000   0.957   2.600  1.00 10.00           H
HETATM    4  O   HOH W   3       0.000   2.600   0.000  1.00 10.00           O
HETATM    5  O   HOH W   4       0.000   2.600   2.600  1.00 10.00           O
HETATM    6  H1  HOH W   4       0.000   3.557   2.600  1.00 10.00           H
HETATM    7  O   HOH W   5       2.600   0.000   0.000  1.00 10.00           O
HETATM    8  O   HOH W   6       2.600   0.000   2.600  1.00 10.00           O
HETATM    9  H1  HOH W   6       2.600   0.957   2.600  1.00 10.00           H
HETATM   10  O   HOH W   7       2.600   2.600   0.000  1.00 10.00           O
HETATM   11  O   HOH W   8       2.600   2.600   2.600  1.00 10.00           O
HETATM   12  H1  HOH W   8       2.600   3.557   2.600  1.00 10.00           H
END
"""


def _hierarchy(pdb_str):
  """Build a hierarchy from an inline PDB string.

  Parameters
  ----------
  pdb_str : str
      PDB-format record text.

  Returns
  -------
  iotbx.pdb.hierarchy.root
      The constructed hierarchy.
  """
  return iotbx.pdb.input(
    source_info=None, lines=pdb_str.split("\n")).construct_hierarchy()


def _water_atoms(hier):
  """Pull the O and placed H of the single HOH in a hierarchy.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Hierarchy containing exactly one HOH residue.

  Returns
  -------
  tuple
      ``(o, hs)``: the O atom and a ``{name: atom}`` dict of its H/D.
  """
  water = [a for a in hier.atoms() if a.parent().resname.strip() == "HOH"]
  o = next(a for a in water if a.element.strip().upper() == "O")
  hs = {a.name.strip(): a for a in water
        if a.element.strip().upper() in ("H", "D")}
  return o, hs


def _unit(a, b):
  """Unit vector from atom/point ``b`` to atom/point ``a``.

  Parameters
  ----------
  a, b : iotbx.pdb.hierarchy.atom or sequence of float
      Endpoints, each either an atom (``.xyz`` is read) or an ``(x, y, z)``.

  Returns
  -------
  scitbx.matrix.col
      The unit vector ``(a - b)``.
  """
  va = matrix.col(a.xyz) if hasattr(a, "xyz") else matrix.col(a)
  vb = matrix.col(b.xyz) if hasattr(b, "xyz") else matrix.col(b)
  return (va - vb).normalize()


def exercise_acceptor_directed():
  """Both protons point at the flanking acceptors, clash-free, with the
  right O-H length and H-O-H angle."""
  hier = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  atoms = list(hier.atoms())
  o, hs = _water_atoms(hier)
  assert set(hs) == {"H1", "H2"}, f"expected H1 and H2; got {sorted(hs)}"

  acc = sorted([a for a in atoms if a.parent().resname.strip() in ("ACA", "ACB")],
               key=lambda a: a.distance(o))
  a1, a2 = acc  # nearer first -> matches H1
  assert _unit(hs["H1"], o).dot(_unit(a1, o)) > 0.9, (
    "H1 should point at the nearest acceptor")
  assert _unit(hs["H2"], o).dot(_unit(a2, o)) > 0.9, (
    "H2 should point at the second acceptor, not open space")

  non_water = [a for a in atoms if a.parent().resname.strip() != "HOH"]
  for h in (hs["H1"], hs["H2"]):
    d = min(h.distance(a) for a in non_water)
    assert d >= wp._WATER_MIN_CLEARANCE - 1e-6, (
      f"placed H too close to a non-water atom: {d:.3f} A")

  oh = (matrix.col(hs["H1"].xyz) - matrix.col(o.xyz)).length()
  assert abs(oh - wp._WATER_OH_XRAY) < 1e-3, f"O-H length off: {oh:.3f}"
  ang = math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))
  assert abs(ang - wp._WATER_HOH_DEG) < 1.0, f"H-O-H angle off: {ang:.1f}"


def exercise_cation_repulsion():
  """Both protons of a metal-coordinating water sit in the hemisphere away
  from the cation, even though the only acceptor is on the metal side."""
  hier = _hierarchy(_MG_WATER_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  mg = next(a for a in hier.atoms() if a.element.strip().upper() == "MG")
  o, hs = _water_atoms(hier)
  assert len(hs) == 2, f"expected two placed H; got {sorted(hs)}"

  to_mg = _unit(mg, o)
  for h in hs.values():
    proj = (matrix.col(h.xyz) - matrix.col(o.xyz)).dot(to_mg)
    assert proj <= 1e-3, (
      f"water H should point away from the cation (proj={proj:.2f})")


def exercise_idempotent():
  """A water that already has two H is left exactly as-is."""
  hier = _hierarchy(_PROTONATED_WATER_PDB)
  before = [(a.name.strip(), tuple(a.xyz)) for a in hier.atoms()]
  wp.place_water_hydrogens(hier)
  after = [(a.name.strip(), tuple(a.xyz)) for a in hier.atoms()]
  assert after == before, (
    f"already-protonated water must be untouched:\n{before}\n{after}")


def exercise_protonated_n_not_acceptor():
  """An N that already carries an H is a donor, not an acceptor: H1 skips
  the nearer protonated N and points at the farther bare N instead."""
  hier = _hierarchy(_DONOR_N_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  o_xyz = matrix.col(o.xyz)
  acc = matrix.col((-1.0, 0.0, 0.0))   # toward the bare N (acceptor)
  don = matrix.col((1.0, 0.0, 0.0))    # toward the protonated N (donor)
  best_acc = max((matrix.col(h.xyz) - o_xyz).normalize().dot(acc)
                 for h in hs.values())
  best_don = max((matrix.col(h.xyz) - o_xyz).normalize().dot(don)
                 for h in hs.values())
  assert best_acc > 0.9, "H should point at the bare (acceptor) N"
  assert best_don < 0.9, "no H should point at the protonated (donor) N"


def exercise_lone_pair_directed():
  """``lone_pair_directed`` aims the O-H at a carbonyl O's sp2 lone-pair
  lobe (~120 deg from C=O); the default aims at the nucleus (~180 deg)."""
  o_c = matrix.col((0.0, 0.0, 0.0))    # carbonyl O
  c = matrix.col((-1.220, 0.0, 0.0))   # its carbon

  def co_h_angle(lone_pair):
    hier = _hierarchy(_CARBONYL_PDB)
    wp.place_water_hydrogens(hier, n_refine=0, lone_pair_directed=lone_pair)
    _, hs = _water_atoms(hier)
    h = min(hs.values(), key=lambda a: (matrix.col(a.xyz) - o_c).length())
    return (c - o_c).angle(matrix.col(h.xyz) - o_c, deg=True)

  off = co_h_angle(False)
  on = co_h_angle(True)
  assert off > 135.0, f"default should aim near the nucleus (got {off:.1f} deg)"
  assert abs(on - 120.0) < 15.0, (
    f"lone-pair placement should give a ~120 deg C=O...H angle (got {on:.1f})")


def exercise_joint_balances_acceptors():
  """With two acceptors too acute (60 deg) for both H-bonds from an
  on-acceptor H1, joint placement balances them -- the worse of the two
  H-bonds is markedly better than greedy's leftover."""
  a = matrix.col((2.800, 0.000, 0.000))
  b = matrix.col((1.400, 2.424, 0.000))

  def worse_alignment(joint):
    hier = _hierarchy(_TWO_ACCEPTOR_60_PDB)
    wp.place_water_hydrogens(hier, n_refine=0, joint=joint)
    o, hs = _water_atoms(hier)
    da = (a - matrix.col(o.xyz)).normalize()
    db = (b - matrix.col(o.xyz)).normalize()
    align_a = max(_unit(h, o).dot(da) for h in hs.values())
    align_b = max(_unit(h, o).dot(db) for h in hs.values())
    return min(align_a, align_b)        # quality of the worse H-bond

  greedy = worse_alignment(False)
  joint = worse_alignment(True)
  assert joint > greedy + 0.1, (
    f"joint should balance the two H-bonds (worse-bond align "
    f"{greedy:.2f} -> {joint:.2f})")
  assert joint > 0.85, (
    f"joint's worse H-bond should still be good (got {joint:.2f})")


def exercise_reorient_existing():
  """``existing_h="reorient"`` strips the existing (mis-oriented) H and
  re-places them toward the acceptor; the default leaves them as-is."""
  acc = matrix.col((2.700, 0.000, 0.000))

  def best_align(hier):
    o, hs = _water_atoms(hier)
    return max((matrix.col(h.xyz) - matrix.col(o.xyz)).normalize().dot(
                 (acc - matrix.col(o.xyz)).normalize()) for h in hs.values())

  kept = _hierarchy(_BAD_PROTONATED_PDB)
  wp.place_water_hydrogens(kept, n_refine=0)
  assert best_align(kept) < 0.5, (
    "default run must not reorient existing H toward the acceptor")

  redone = _hierarchy(_BAD_PROTONATED_PDB)
  wp.place_water_hydrogens(redone, n_refine=0, existing_h="reorient")
  _, hs = _water_atoms(redone)
  assert len(hs) == 2, f"reorient must leave exactly two H; got {sorted(hs)}"
  assert best_align(redone) > 0.9, (
    'existing_h="reorient" should re-place an H toward the acceptor')


def exercise_single_h_water():
  """A water with one H is a possible hydroxide: untouched and reported by
  default, and completed on the *deposited* proton's cone on request."""
  acc = matrix.col((2.700, 0.000, 0.000))
  h1_in = matrix.col((0.000, 0.957, 0.000))

  kept = _hierarchy(_SINGLE_H_WATER_PDB)
  res = wp.place_water_hydrogens(kept, n_refine=0)
  o, hs = _water_atoms(kept)
  assert len(hs) == 1, (
    f"default must not complete a single-H water; got {sorted(hs)}")
  assert (matrix.col(hs["H1"].xyz) - h1_in).length() < 1e-6, (
    "the deposited proton must not move")
  assert [r[2] for r in res.partial_waters] == ["kept"]

  done = _hierarchy(_SINGLE_H_WATER_PDB)
  res = wp.place_water_hydrogens(done, n_refine=0, existing_h="complete")
  o, hs = _water_atoms(done)
  assert len(hs) == 2, f"complete must add the partner; got {sorted(hs)}"
  assert (matrix.col(hs["H1"].xyz) - h1_in).length() < 1e-6, (
    "completing must not move the deposited proton")
  O = matrix.col(o.xyz)
  ang = math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))
  assert abs(ang - wp._WATER_HOH_DEG) < 1e-3, (
    f"H-O-H must be canonical against the deposited H; got {ang:.2f} deg")
  assert abs((matrix.col(hs["H2"].xyz) - O).length()
             - wp._WATER_OH_XRAY) < 1e-6, "new O-H length"
  assert _unit(hs["H2"], o).dot((acc - O).normalize()) > 0.9, (
    "the new proton should aim at the acceptor")
  assert [r[2] for r in res.partial_waters] == ["completed"]

  redone = _hierarchy(_SINGLE_H_WATER_PDB)
  res = wp.place_water_hydrogens(redone, n_refine=0, existing_h="reorient")
  o, hs = _water_atoms(redone)
  assert len(hs) == 2, f"reorient must give two H; got {sorted(hs)}"
  assert (matrix.col(hs["H1"].xyz) - h1_in).length() > 0.1, (
    "reorient must re-place the deposited proton")
  assert [r[2] for r in res.partial_waters] == ["stripped"]


def exercise_single_h_metal_annotation():
  """The single-H report flags a coordinating cation, using a first-shell
  cutoff rather than the looser proton-repulsion radius."""
  near = wp.place_water_hydrogens(
    _hierarchy(_MG_SINGLE_H_PDB), n_refine=0).partial_waters
  assert len(near) == 1, f"expected one single-H water; got {near}"
  rid, metal, action = near[0]
  assert rid == "HOH W 1", rid
  assert action == "kept"
  assert metal is not None, "a 2.1 A Mg should be reported as coordinating"
  assert metal[0] == "MG" and abs(metal[1] - 2.100) < 1e-3, metal

  far = wp.place_water_hydrogens(
    _hierarchy(_MG_SINGLE_H_FAR_PDB), n_refine=0).partial_waters
  assert len(far) == 1 and far[0][1] is None, (
    f"a 2.9 A Mg is inside the repulsion radius but is not a bond; got {far}")


def exercise_hetatm_flag():
  """Placed H inherit the parent O's HETATM flag, so a HETATM water does
  not acquire ATOM-record protons."""
  hier = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  assert o.hetero, "fixture water should be a HETATM"
  for name, a in sorted(hs.items()):
    assert a.hetero == o.hetero, (
      f"{name} must match the parent O record type")


def exercise_heavy_cation_repulsion():
  """The cation set is not limited to the first row: a Pt keeps both
  protons out of its hemisphere just as an Mg does."""
  hier = _hierarchy(_MG_WATER_PDB.replace("MG", "PT"))
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  O = matrix.col(o.xyz)
  cd = (matrix.col((2.500, 0.000, 0.000)) - O).normalize()
  for name, h in sorted(hs.items()):
    assert (matrix.col(h.xyz) - O).dot(cd) <= 1e-9, (
      f"{name} must stay out of the Pt hemisphere")


def exercise_refinement_reduces_clashes():
  """Refinement relaxes the water-water H clashes the greedy pass leaves in
  a tight cluster."""
  greedy = _hierarchy(_WATER_CLUSTER_PDB)
  wp.place_water_hydrogens(greedy, n_refine=0, oh_length=wp._WATER_OH_NEUTRON)
  refined = _hierarchy(_WATER_CLUSTER_PDB)
  wp.place_water_hydrogens(refined, n_refine=5, oh_length=wp._WATER_OH_NEUTRON)

  n_greedy = wp._water_clash_stats(greedy)[1]
  n_refined = wp._water_clash_stats(refined)[1]
  assert n_greedy >= 1, (
    f"cluster should clash without refinement (got {n_greedy})")
  assert n_refined < n_greedy, (
    f"refinement should reduce close contacts ({n_greedy} -> {n_refined})")


def exercise_completed_water_survives_refinement():
  """Refinement sweeps and basin kicks must hold a completed water's O-H1
  axis fixed: the deposited proton never moves and the partner stays on its
  cone, instead of the axis being re-derived from the environment."""
  hier = _hierarchy(_MIXED_CLUSTER_PDB)
  before = {}
  for ag in hier.atom_groups():
    for a in ag.atoms():
      if a.element.strip().upper() == "H":
        before[(ag.parent().resseq.strip(), a.name.strip())] = tuple(a.xyz)
  assert len(before) == 4, f"fixture should deposit four H; got {len(before)}"

  states = []
  res = wp.place_water_hydrogens(
    hier, n_refine=3, n_basin=2, existing_h="complete",
    oh_length=wp._WATER_OH_NEUTRON,
    on_state=lambda label, stats: states.append(stats[1]))
  assert len(res.partial_waters) == 4
  # The point of the fixture: refinement and basin-hopping really do run.
  assert len(states) > 1 and max(states) > 0, (
    f"fixture must clash so the sweeps engage; states={states}")

  for ag in hier.atom_groups():
    o = next(a for a in ag.atoms() if a.element.strip().upper() == "O")
    hs = [a for a in ag.atoms() if a.element.strip().upper() == "H"]
    assert len(hs) == 2, f"every water should end with two H; got {len(hs)}"
    for a in hs:
      key = (ag.parent().resseq.strip(), a.name.strip())
      if key in before:
        moved = (matrix.col(a.xyz) - matrix.col(before[key])).length()
        assert moved < 1e-9, f"deposited {key} moved {moved:.3f} A"
    ang = math.degrees(_unit(hs[0], o).angle(_unit(hs[1], o)))
    assert abs(ang - wp._WATER_HOH_DEG) < 0.1, (
      f"H-O-H {ang:.2f} deg on water {ag.parent().resseq.strip()}")


def exercise_element_override():
  """``element="D"`` forces deuterium (named D1/D2); the default is H."""
  default = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(default, n_refine=0)
  _, hs = _water_atoms(default)
  assert set(hs) == {"H1", "H2"}, f"default should place H; got {sorted(hs)}"

  deut = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(deut, n_refine=0, element="D")
  _, hs = _water_atoms(deut)
  assert set(hs) == {"D1", "D2"}, f"element='D' should place D; got {sorted(hs)}"
  for d in hs.values():
    assert d.element.strip().upper() == "D", "placed atom element must be D"


def exercise_oh_length_auto():
  """With no ``oh_length`` the placer picks the canonical distance from the
  model: neutron where D is present, X-ray otherwise. An explicit value
  overrides the choice in either direction."""

  def oh(hier):
    o, hs = _water_atoms(hier)
    lengths = [(matrix.col(h.xyz) - matrix.col(o.xyz)).length()
               for h in hs.values()]
    assert max(lengths) - min(lengths) < 1e-6, lengths
    return lengths[0]

  # Hydrogenous model: X-ray.
  xray = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(xray, n_refine=0)
  assert abs(oh(xray) - wp._WATER_OH_XRAY) < 1e-6, (
    f"an H-only model should get the X-ray length; got {oh(xray):.3f}")

  # A D anywhere in the model means neutron, even on another residue.
  deut = _hierarchy(_TWO_ACCEPTOR_PDB.replace(
    "HETATM    3  O   ACB D   2       4.299   5.000   7.711  1.00 10.00           O",
    "HETATM    3  D   ACB D   2       4.299   5.000   7.711  1.00 10.00           D"))
  wp.place_water_hydrogens(deut, n_refine=0)
  assert abs(oh(deut) - wp._WATER_OH_NEUTRON) < 1e-6, (
    f"a model carrying D should get the neutron length; got {oh(deut):.3f}")

  # An explicit value wins over the heuristic.
  forced = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(forced, n_refine=0,
                           oh_length=wp._WATER_OH_NEUTRON)
  assert abs(oh(forced) - wp._WATER_OH_NEUTRON) < 1e-6, (
    "an explicit oh_length must override the heuristic")


def exercise_detect_neutron():
  """``_detect_neutron`` prefers the experiment type, then D-atom presence."""
  def pdb_in(lines):
    return iotbx.pdb.input(source_info=None, lines=lines.split("\n"))

  neutron = "EXPDTA    NEUTRON DIFFRACTION\n" + _TWO_ACCEPTOR_PDB
  xray = "EXPDTA    X-RAY DIFFRACTION\n" + _TWO_ACCEPTOR_PDB
  # No experiment record, but a D atom present -> neutron by fallback.
  d_atom = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  D   DNR E   1       3.600   0.000   0.000  1.00 10.00           D
END
"""

  for src, want in ((neutron, True), (xray, False), (d_atom, True)):
    pi = pdb_in(src)
    got, _ = wp._detect_neutron(pi, pi.construct_hierarchy())
    assert got is want, f"_detect_neutron returned {got}, expected {want}"


def run():
  """Run every exercise and print the CPU times and ``OK`` on success."""
  exercise_acceptor_directed()
  exercise_cation_repulsion()
  exercise_protonated_n_not_acceptor()
  exercise_lone_pair_directed()
  exercise_joint_balances_acceptors()
  exercise_idempotent()
  exercise_reorient_existing()
  exercise_single_h_water()
  exercise_single_h_metal_annotation()
  exercise_hetatm_flag()
  exercise_heavy_cation_repulsion()
  exercise_refinement_reduces_clashes()
  exercise_completed_water_survives_refinement()
  exercise_element_override()
  exercise_oh_length_auto()
  exercise_detect_neutron()
  print(format_cpu_times())
  print("OK")


if __name__ == "__main__":
  run()
