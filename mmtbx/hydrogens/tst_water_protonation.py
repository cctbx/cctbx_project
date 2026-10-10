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

* **H2 aims at a reachable acceptor.** H2 is confined to the H-O-H cone about
  O-H1, so an acceptor well off that cone cannot be donated to. The nearest
  free acceptor being unreachable must not stop H2 from aiming at a further
  one that lies on the cone.

* **Refinement.** A tight cluster of bare waters clashes after the greedy
  pass; the relaxation sweeps must reduce the count of close H-H contacts.

* **A completed water stays completed.** In a clashing cluster where half
  the waters carry a deposited proton, the refinement sweeps and basin
  kicks must leave every deposited proton where it is and keep each new
  partner on its cone, rather than re-deriving the O-H1 axis.

* **Element override.** ``element="D"`` forces deuterium (named ``D1``/``D2``)
  onto an HOH water; the default places H.

* **Environment hydrogen count.** ``count_environment_hydrogens`` counts H and
  D outside water residues and ignores the waters themselves, so a model whose
  only hydrogens sit on its waters reads as zero.

* **Experiment detection.** ``_detect_neutron`` reads the experiment type
  (EXPDTA / ``_exptl.method``) and falls back to the presence of D atoms.

* **O-H length auto-selection.** With no ``oh_length`` the placer takes the
  neutron distance for a model carrying D and the X-ray distance otherwise;
  an explicit value overrides both.

* **Crystal symmetry.** An acceptor reachable only across a cell face draws
  an H once a crystal symmetry is given, and is invisible without one; the
  symmetry-related copies leave the model's own atoms where they are.

* **Protons across a lattice contact.** Two waters donating across a lattice
  contact keep their protons clear of each other's symmetry-related protons,
  not just of the oxygens they aim at.

* **Electron microscopy.** The program treats an electron microscopy model as
  isolated: its cell is the map box, not a lattice.

* **Reorient keeps the isotope.** ``existing_h="reorient"`` with the automatic
  element re-places a water's protons as the element it carried (D on an
  HOH, H on a DOD), and the program warns when that puts D at the X-ray
  O-H length.

* **Missing element columns.** A model with atoms lacking an element symbol
  is refused, by the placer and the program, rather than having every water
  skipped.

* **Multi-model input.** The program and the placer reject more than one
  MODEL, whose copies would otherwise share one environment.

* **Output file name.** The default is ``<model-stem>_waters_protonated``;
  ``output.prefix``, ``output.suffix`` and ``output.serial`` change it.
"""

from __future__ import absolute_import, division, print_function

import math
import os
from io import StringIO

import iotbx.cif
import iotbx.pdb
import libtbx.load_env
from iotbx.cli_parser import CCTBXParser, run_program
from iotbx.data_manager import DataManager
from scitbx import matrix
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import Sorry, format_cpu_times, null_out

from mmtbx.hydrogens import water_protonation as wp
from mmtbx.programs import water_protonation as wp_program


# One HOH (bare O) with two acceptor O placed on the 104.5 deg H-O-H cone:
# A1 at 2.6 A along +x, A2 at 2.8 A at 104.5 deg from +x (+z azimuth).  The
# placer should orient H1 -> A1 (nearer) and H2 -> A2.
_TWO_ACCEPTOR_PDB = """\
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  O   ACA D   1       7.600   5.000   5.000  1.00 10.00           O
HETATM    3  O   ACB D   2       4.299   5.000   7.711  1.00 10.00           O
END
"""

# The same, without element columns, as older PDB files are written.
_NO_ELEMENT_PDB = "\n".join(l[:76].rstrip()
                            for l in _TWO_ACCEPTOR_PDB.split("\n"))

# A water by an unidentified atom in a crystal. The wwPDB writes such an
# atom as UNX with element X, which has no scattering type.
_UNKNOWN_ATOM_PDB = """\
CRYST1   20.000   20.000   20.000  90.00  90.00  90.00 P 1
HETATM    1  UNK UNX A   1       5.000   5.000   5.000  1.00 10.00           X
HETATM    2  O   HOH W   1       7.700   5.000   5.000  1.00 10.00           O
END
"""

# Water O coordinating an Mg (2.5 A along +x). The ONLY acceptor sits on
# the Mg side (30 deg off the O->Mg axis), so acceptor-direction alone
# would pull an H toward the metal; only the cation repulsion keeps both
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

# A water that is already protonated (O + 2 H) must be left untouched.
_PROTONATED_WATER_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       0.000   0.957   0.000  1.00 10.00           H
HETATM    3  H2  HOH W   1       0.926  -0.239   0.000  1.00 10.00           H
END
"""

# A protonated water whose H point AWAY from the lone acceptor (+x): with
# existing_h="reorient" the H are stripped and re-placed toward the acceptor.
_BAD_PROTONATED_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  H1  HOH W   1      -0.984   0.000   0.000  1.00 10.00           H
HETATM    3  H2  HOH W   1       0.000  -0.984   0.000  1.00 10.00           H
HETATM    4  O   ACA D   1       2.700   0.000   0.000  1.00 10.00           O
END
"""

# The same mis-oriented water carrying D instead of H.
_BAD_DEUTERATED_PDB = (_BAD_PROTONATED_PDB
                       .replace(" H1 ", " D1 ").replace(" H2 ", " D2 ")
                       .replace("           H\n", "           D\n"))

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

# The same water by the -a face of a 10 A cell, the Mg across it: only the
# Mg's -a translate, 2.1 A from the O, coordinates the water.
_MG_SINGLE_H_SYM_PDB = """\
CRYST1   10.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1 MG    MG A   1       8.900  15.000  15.000  1.00 10.00          MG
HETATM    2  O   HOH W   1       1.000  15.000  15.000  1.00 10.00           O
HETATM    3  H1  HOH W   1       1.000  15.957  15.000  1.00 10.00           H
END
"""

# H1 goes to ACA (+x, nearest). ACB is the next nearest but sits only 30 deg
# off the O-H1 axis, so no point on the 104.5 deg cone can aim at it. ACC is
# further away yet lies exactly on the cone, so it is the one H2 can donate to.
_UNREACHABLE_ACCEPTOR_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   ACA D   1       2.700   0.000   0.000  1.00 10.00           O
HETATM    3  O   ACB E   1       2.425   1.400   0.000  1.00 10.00           O
HETATM    4  O   ACC F   1      -0.725  -2.807   0.000  1.00 10.00           O
END
"""

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


# A protonated residue (one H, one D) beside a water that also carries H.
# Only the residue's two count.
_ENV_H_PDB = """\
ATOM      1  N   ASN A   1      -1.500   0.000   0.000  1.00 20.00           N
ATOM      2  H   ASN A   1      -2.000   0.800   0.000  1.00 20.00           H
ATOM      3  CG  ASN A   1      -2.500   0.500   0.000  1.00 20.00           C
ATOM      4  D   ASN A   1      -3.000   1.200   0.000  1.00 20.00           D
HETATM    5  O   HOH W   1       1.500   0.000   0.000  1.00 25.00           O
HETATM    6  H1  HOH W   1       0.516   0.000   0.000  1.00 25.00           H
HETATM    7  H2  HOH W   1       1.746   0.953   0.000  1.00 25.00           H
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


def _hierarchy_and_symmetry(pdb_str):
  """Build a hierarchy and its crystal symmetry from an inline PDB string.

  Parameters
  ----------
  pdb_str : str
      PDB-format record text, CRYST1 included.

  Returns
  -------
  tuple
      ``(hierarchy, crystal_symmetry)``.
  """
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  return inp.construct_hierarchy(), inp.crystal_symmetry()


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


def _fallback_cage_pdb():
  """A bare water whose every H1 candidate clashes.

  A C atom sits 2.42 A from the water O along each of the placer's fallback
  directions, so every fallback candidate for H1 lies 1.463 A from one. The
  one acceptor, an O 2.8 A out, lies 6 deg off fallback direction 4; its
  candidate clashes too, but at 1.472 A, so the fallback settles on it. The
  gaps between the C atoms leave parts of the H2 cone clear.

  Returns
  -------
  str
      PDB-format record text.
  """
  o = matrix.col((10.0, 10.0, 10.0))
  dirs = [matrix.col(d).normalize() for d in wp._WATER_FALLBACK_DIRECTIONS]
  acc = dirs[4].rotate_around_origin(axis=dirs[4].ortho(), angle=6.0, deg=True)
  sites = [("O   HOH W", 1, o), ("O   ACA D", 1, o + acc * 2.8)]
  sites += [("C   BLK X", i + 1, o + d * 2.42) for i, d in enumerate(dirs)]
  return "".join(
    f"HETATM{k + 1:5d}  {label}{resseq:4d}    {x:8.3f}{y:8.3f}{z:8.3f}"
    f"  1.00 10.00           {label[0]}\n"
    for k, (label, resseq, (x, y, z)) in enumerate(sites)) + "END\n"


# Water 1.3 A from the +a face of a 10 A cell, its only acceptor 7.2 A away
# across the cell. The acceptor's +a translate sits 2.8 A from the water, so
# the H-bond exists only once crystal symmetry is honoured.
_LATTICE_CONTACT_PDB = """\
CRYST1   10.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH W   1       8.700  15.000  15.000  1.00 10.00           O
HETATM    2  O   ACA D   1       1.500  15.000  15.000  1.00 10.00           O
END
"""


# Two waters 7.2 A apart in a 10 A cell, so each one's only acceptor in range
# is the other's lattice translate, 2.8 A away on the far side. Protons aimed
# straight at those translates collide with each other at 0.89 A.
_CONTACT_PAIR_PDB = """\
CRYST1   10.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH W   1       1.000  15.000  15.000  1.00 10.00           O
HETATM    2  O   HOH W   2       8.200  15.000  15.000  1.00 10.00           O
END
"""


# Three waters that already carry H, in P-1. No H-H contact inside the
# model; across symmetry, H1 of W 2 lies 1.486 A from H1 of the -a translate
# of W 3, and H1 of W 1 lies 1.600 A from its own equivalent through the
# inversion centre at (0, 15, 15).
_SYM_CONTACT_PDB = """\
CRYST1   10.000   30.000   30.000  90.00  90.00  90.00 P -1
HETATM    1  O   HOH W   1       1.757  15.000  15.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       0.800  15.000  15.000  1.00 10.00           H
HETATM    3  H2  HOH W   1       1.997  15.927  15.000  1.00 10.00           H
HETATM    4  O   HOH W   2       1.000   7.500   7.500  1.00 10.00           O
HETATM    5  H1  HOH W   2       0.043   7.500   7.500  1.00 10.00           H
HETATM    6  H2  HOH W   2       1.240   8.427   7.500  1.00 10.00           H
HETATM    7  O   HOH W   3       7.600   7.500   7.500  1.00 10.00           O
HETATM    8  H1  HOH W   3       8.557   7.500   7.500  1.00 10.00           H
HETATM    9  H2  HOH W   3       7.360   8.427   7.500  1.00 10.00           H
END
"""


# A water on the two-fold axis of P2, carrying one H. The two-fold maps the
# O onto itself and the H onto the water's second H site.
_ON_TWO_FOLD_PDB = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  O   HOH W   1       0.000   5.000   0.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       0.757   5.586   0.000  1.00 10.00           H
END
"""


# A water split between altlocs A and B, 1.39 A apart, and a blank acceptor
# 2.8 A from conformer A along +x. Conformer B's O sits 1.23 A from the H
# conformer A would aim at the acceptor, but the two conformers never
# coexist.
_SPLIT_WATER_PDB = """\
HETATM    1  O  AHOH W   1       5.000   5.000   5.000  0.50 10.00           O
HETATM    2  O  BHOH W   1       5.700   6.200   5.000  0.50 10.00           O
HETATM    3  O   ACA D   1       7.800   5.000   5.000  1.00 10.00           O
END
"""

# One H name both blank and in altloc A: an improper altloc.
_IMPROPER_ALTLOC_PDB = """\
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       5.957   5.000   5.000  1.00 10.00           H
HETATM    3  H1 AHOH W   1       5.000   5.957   5.000  0.50 10.00           H
END
"""


# A blank O whose D sit in altlocs A (0.6) and B (0.4), two in each.
_ALTLOC_D_WATER_PDB = """\
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  D1 AHOH W   1       5.984   5.000   5.000  0.60 10.00           D
HETATM    3  D2 AHOH W   1       4.754   5.953   5.000  0.60 10.00           D
HETATM    4  D1 BHOH W   1       5.000   5.000   5.984  0.40 10.00           D
HETATM    5  D2 BHOH W   1       5.953   5.000   4.754  0.40 10.00           D
END
"""

# An O split between altlocs A and B 0.3 A apart, sharing one blank D.
_ALTLOC_O_WATER_PDB = """\
HETATM    1  O  AHOH W   1       5.000   5.000   5.000  0.50 10.00           O
HETATM    2  O  BHOH W   1       5.300   5.000   5.000  0.50 10.00           O
HETATM    3  D1  HOH W   1       5.150   5.970   5.000  1.00 10.00           D
END
"""

# One residue holding a water in altloc A and a sulfate in altloc B.
_WATER_SO4_ALTLOC_PDB = """\
HETATM    1  O  AHOH W   1       5.000   5.000   5.000  0.50 10.00           O
HETATM    2  S  BSO4 W   1       7.000   5.000   5.000  0.50 10.00           S
HETATM    3  O1 BSO4 W   1       7.000   6.430   5.000  0.50 10.00           O
END
"""

# W 2 blank and W 3 split between altlocs A and B 0.6 A apart, all carrying
# H, in P1: H1 of W 2 lies 1.486 A from the -a translate of W 3's H1 in A
# and 1.603 A from that of its H1 in B.
_SPLIT_SYM_CONTACT_PDB = """\
CRYST1   10.000   30.000   30.000  90.00  90.00  90.00 P 1
HETATM    1  O   HOH W   2       1.000   7.500   7.500  1.00 10.00           O
HETATM    2  H1  HOH W   2       0.043   7.500   7.500  1.00 10.00           H
HETATM    3  H2  HOH W   2       1.240   8.427   7.500  1.00 10.00           H
HETATM    4  O  AHOH W   3       7.600   7.500   7.500  0.50 10.00           O
HETATM    5  H1 AHOH W   3       8.557   7.500   7.500  0.50 10.00           H
HETATM    6  H2 AHOH W   3       7.360   8.427   7.500  0.50 10.00           H
HETATM    7  O  BHOH W   3       7.600   7.500   8.100  0.50 10.00           O
HETATM    8  H1 BHOH W   3       8.557   7.500   8.100  0.50 10.00           H
HETATM    9  H2 BHOH W   3       7.360   6.573   8.100  0.50 10.00           H
END
"""


# The same water and acceptor in two models.
_MULTI_MODEL_PDB = """\
MODEL        1
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  O   ACA D   1       7.600   5.000   5.000  1.00 10.00           O
ENDMDL
MODEL        2
HETATM    1  O   HOH W   1       5.010   5.000   5.000  1.00 10.00           O
HETATM    2  O   ACA D   1       7.600   5.000   5.000  1.00 10.00           O
ENDMDL
END
"""


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

  oh_target, hoh_target = wp._WATER_GEOMETRY["xray"]
  oh = (matrix.col(hs["H1"].xyz) - matrix.col(o.xyz)).length()
  assert abs(oh - oh_target) < 1e-3, f"O-H length off: {oh:.3f}"
  ang = math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))
  assert abs(ang - hoh_target) < 1.0, f"H-O-H angle off: {ang:.1f}"


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
  """A water that already has two H is left exactly as-is, as is a blank O
  carrying two D in each of altlocs A and B."""
  for pdb_str in (_PROTONATED_WATER_PDB, _ALTLOC_D_WATER_PDB):
    hier = _hierarchy(pdb_str)
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
  lobe (~120 deg from C=O); the default aims at the nucleus (~180 deg).
  With the carbonyl O split between altlocs 0.25 A apart, each copy's lobes
  still come from its C alone, not from a bond to the other copy."""
  c = matrix.col((-1.220, 0.0, 0.0))   # the carbonyl C
  split = _CARBONYL_PDB.replace(
    "HETATM    2  O   ACO A   1       0.000   0.000   0.000  1.00 10.00",
    "HETATM    2  O  AACO A   1       0.000   0.000   0.000  0.50 10.00").replace(
    "END", "HETATM    6  O  BACO A   1       0.250   0.000   0.000  0.50 10.00"
    "           O\nEND")

  def co_h_angle(pdb_str, lone_pair):
    hier = _hierarchy(pdb_str)
    wp.place_water_hydrogens(hier, n_refine=0, lone_pair_directed=lone_pair)
    _, hs = _water_atoms(hier)
    acc = [matrix.col(a.xyz) for a in hier.atoms()
           if a.parent().resname == "ACO" and a.element.strip() == "O"]
    h, o_c = min(((matrix.col(a.xyz), o) for a in hs.values() for o in acc),
                 key=lambda t: (t[0] - t[1]).length())
    return (c - o_c).angle(h - o_c, deg=True)

  off = co_h_angle(_CARBONYL_PDB, False)
  assert off > 135.0, f"default should aim near the nucleus (got {off:.1f} deg)"
  for pdb_str in (_CARBONYL_PDB, split):
    on = co_h_angle(pdb_str, True)
    assert abs(on - 120.0) < 15.0, (
      f"lone-pair placement should give a ~120 deg C=O...H angle "
      f"(got {on:.1f})")


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
  oh_target, hoh_target = wp._WATER_GEOMETRY["xray"]
  ang = math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))
  assert abs(ang - hoh_target) < 1e-3, (
    f"H-O-H must be canonical against the deposited H; got {ang:.2f} deg")
  assert abs((matrix.col(hs["H2"].xyz) - O).length()
             - oh_target) < 1e-6, "new O-H length"
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
  cutoff rather than the looser proton-repulsion radius, and a cation that
  only crystal symmetry brings next to the water even when no water is
  placed."""
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

  hier, cs = _hierarchy_and_symmetry(_MG_SINGLE_H_SYM_PDB)
  sym = wp.place_water_hydrogens(
    hier, n_refine=0, crystal_symmetry=cs).partial_waters
  assert sym[0][1] is not None and abs(sym[0][1][1] - 2.100) < 1e-3, sym


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


def exercise_h2_reachable_acceptor():
  """H2 aims at an acceptor it can actually reach on the H-O-H cone, not at
  the nearest free one when that lies off the cone."""
  hier = _hierarchy(_UNREACHABLE_ACCEPTOR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  assert len(hs) == 2, sorted(hs)
  O = matrix.col(o.xyz)
  aca = (matrix.col((2.700, 0.000, 0.000)) - O).normalize()
  acc = (matrix.col((-0.725, -2.807, 0.000)) - O).normalize()
  d1, d2 = _unit(hs["H1"], o), _unit(hs["H2"], o)
  assert d1.dot(aca) > 0.99, f"H1 should take the nearest acceptor; got {d1.dot(aca):.3f}"
  assert d2.dot(acc) > 0.9, (
    f"H2 should aim at the reachable acceptor on the cone; "
    f"alignment {d2.dot(acc):.3f}")


def exercise_h2_after_acceptor_fallback():
  """H2 ignores H1's acceptor when the H1 fallback chose it.

  Places on ``_fallback_cage_pdb`` with the gas-phase geometry it is
  built for, where every H1 candidate clashes and the fallback settles on
  the lone acceptor. With no other acceptor, H2 must take the clash-free
  point of the sampled H-O-H cone with the most clearance; scoring H1's
  acceptor, which every cone point sees at the same angle, would leave the
  choice to rounding noise.
  """
  hier = _hierarchy(_fallback_cage_pdb())
  wp.place_water_hydrogens(hier, n_refine=0, geometry="gas_phase")
  o, hs = _water_atoms(hier)
  acc = next(a for a in hier.atoms() if a.parent().resname == "ACA")
  d1 = _unit(hs["H1"], o)
  assert d1.dot(_unit(acc, o)) > 0.9999, "H1 is not on the acceptor"
  heavy = [matrix.col(a.xyz) for a in hier.atoms()
           if a.parent().resname != "HOH"]
  def clearance(pt):
    return min((pt - x).length() for x in heavy)
  # H1 clashes, so it came from the fallback.
  assert clearance(matrix.col(hs["H1"].xyz)) < wp._WATER_MIN_CLEARANCE
  p, q = (matrix.col(v) for v in wp._ortho_frame(d1.elems))
  oh_gas, hoh_gas = wp._WATER_GEOMETRY["gas_phase"]
  hoh = math.radians(hoh_gas)
  n = wp._WATER_CONE_SAMPLES
  cone = [matrix.col(o.xyz) + oh_gas * (
            d1 * math.cos(hoh) + (p * math.cos(2 * math.pi * k / n)
                                  + q * math.sin(2 * math.pi * k / n))
            * math.sin(hoh)) for k in range(n)]
  best = max((c for c in cone if clearance(c) >= wp._WATER_MIN_CLEARANCE),
             key=clearance)
  h2 = matrix.col(hs["H2"].xyz)
  assert (h2 - best).length() < 1e-6, (
    f"H2 clearance {clearance(h2):.3f} A, clearest on the cone "
    f"{clearance(best):.3f} A")


def exercise_refinement_reduces_clashes():
  """Refinement relaxes the water-water H clashes the greedy pass leaves in
  a tight cluster (placed at the 0.984 A O-H and gas-phase H-O-H it is
  built for)."""
  geometry = dict(oh_length=0.984, geometry="gas_phase")
  greedy = _hierarchy(_WATER_CLUSTER_PDB)
  wp.place_water_hydrogens(greedy, n_refine=0, **geometry)
  refined = _hierarchy(_WATER_CLUSTER_PDB)
  wp.place_water_hydrogens(refined, n_refine=5, **geometry)

  n_greedy = wp._water_clash_stats(greedy)[1]
  n_refined = wp._water_clash_stats(refined)[1]
  assert n_greedy >= 1, (
    f"cluster should clash without refinement (got {n_greedy})")
  assert n_refined < n_greedy, (
    f"refinement should reduce close contacts ({n_greedy} -> {n_refined})")
  # The program's contact listing reports the same contacts, closest first.
  listed = [c[0] for c in wp._worst_water_clashes(greedy)]
  assert len(listed) == n_greedy, (len(listed), n_greedy)
  assert listed == sorted(listed), listed


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
    geometry="neutron",
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
    assert abs(ang - wp._WATER_GEOMETRY["neutron"][1]) < 0.1, (
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
  """With no ``geometry`` the placer picks it from the model: neutron where
  D is present, X-ray otherwise, both cctbx's restraint targets for HOH,
  compared with the restraint library when it is installed. An explicit
  ``geometry`` or ``oh_length`` overrides the choice."""

  def oh(hier):
    o, hs = _water_atoms(hier)
    lengths = [(matrix.col(h.xyz) - matrix.col(o.xyz)).length()
               for h in hs.values()]
    assert max(lengths) - min(lengths) < 1e-6, lengths
    return lengths[0]

  def hoh(hier):
    o, hs = _water_atoms(hier)
    return math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))

  xray_oh, xray_hoh = wp._WATER_GEOMETRY["xray"]
  neutron_oh, neutron_hoh = wp._WATER_GEOMETRY["neutron"]
  path = libtbx.env.find_in_repositories(
    relative_path="chem_data/geostd/h/data_HOH.cif", test=os.path.isfile)
  if path is not None:
    hoh_entry = iotbx.cif.reader(file_path=path).model()["comp_HOH"]
    for key, value in (("_chem_comp_bond.value_dist", xray_oh),
                       ("_chem_comp_bond.value_dist_neutron", neutron_oh),
                       ("_chem_comp_angle.value_angle", xray_hoh),
                       ("_chem_comp_angle.value_angle", neutron_hoh)):
      assert [float(v) for v in hoh_entry[key]] == [value] * len(
        hoh_entry[key]), (key, list(hoh_entry[key]), value)

  # Hydrogenous model: X-ray.
  xray = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(xray, n_refine=0)
  assert abs(oh(xray) - xray_oh) < 1e-6, (
    f"an H-only model should get the X-ray length; got {oh(xray):.3f}")
  assert abs(hoh(xray) - xray_hoh) < 1e-6, hoh(xray)

  # A D anywhere in the model means neutron, even on another residue.
  deut_pdb = _TWO_ACCEPTOR_PDB.replace(
    "HETATM    3  O   ACB D   2       4.299   5.000   7.711  1.00 10.00           O",
    "HETATM    3  D   ACB D   2       4.299   5.000   7.711  1.00 10.00           D")
  deut = _hierarchy(deut_pdb)
  wp.place_water_hydrogens(deut, n_refine=0)
  assert abs(oh(deut) - neutron_oh) < 1e-6, (
    f"a model carrying D should get the neutron length; got {oh(deut):.3f}")

  # An explicit value wins over the heuristic.
  forced = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(forced, n_refine=0, oh_length=neutron_oh)
  assert abs(oh(forced) - neutron_oh) < 1e-6, (
    "an explicit oh_length must override the heuristic")

  # So does an explicit geometry.
  gas = _hierarchy(deut_pdb)
  wp.place_water_hydrogens(gas, n_refine=0, geometry="gas_phase")
  assert max(abs(oh(gas) - wp._WATER_GEOMETRY["gas_phase"][0]),
             abs(hoh(gas) - wp._WATER_GEOMETRY["gas_phase"][1])) < 1e-6, (
    oh(gas), hoh(gas))


def exercise_environment_hydrogen_count():
  """``count_environment_hydrogens`` counts H and D outside the waters, and
  is not fooled by a model whose only hydrogens are on its waters."""
  assert wp.count_environment_hydrogens(_hierarchy(_ENV_H_PDB)) == 2

  # Waters keep their H; the residue loses its two.
  bare = "\n".join(l for l in _ENV_H_PDB.split("\n")
                    if " H   ASN" not in l and " D   ASN" not in l)
  hier = _hierarchy(bare)
  assert wp.count_environment_hydrogens(hier) == 0, (
    "hydrogens on waters must not count as environment hydrogens")
  assert sum(1 for a in hier.atoms() if a.element_is_hydrogen()) == 2, (
    "the water H should still be present")

  # DOD is a water alias too, so deuterium on it must not count either.
  dod = "\n".join(
    (l.replace("HOH", "DOD").replace(" H1 ", " D1 ").replace(" H2 ", " D2 ")
      .replace("           H", "           D")) if "HOH" in l else l
    for l in _ENV_H_PDB.split("\n"))
  hier = _hierarchy(dod)
  water_d = [a for a in hier.atoms()
             if a.parent().resname.strip() == "DOD" and a.element.strip() == "D"]
  assert len(water_d) == 2, f"fixture should carry two water D; got {len(water_d)}"
  assert wp.count_environment_hydrogens(hier) == 2

  # Solvent-only model.
  assert wp.count_environment_hydrogens(_hierarchy(_WATER_CLUSTER_PDB)) == 0


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


def exercise_crystal_symmetry():
  """An acceptor reachable only as a lattice translate draws an H once a
  crystal symmetry is given, and is invisible without one."""
  hier, cs = _hierarchy_and_symmetry(_LATTICE_CONTACT_PDB)
  wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs)
  o, hs = _water_atoms(hier)
  aimed = max(_unit(h, o)[0] for h in hs.values())
  assert aimed > 0.99, (
    f"no H aims across the cell face at the translated acceptor "
    f"(best {aimed:.3f})")

  hier, _ = _hierarchy_and_symmetry(_LATTICE_CONTACT_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  isolated = max(_unit(h, o)[0] for h in hs.values())
  assert isolated < 0.9, (
    f"an isolated water should not find the acceptor (best {isolated:.3f})")

  # A translate keeps its altloc: a water in A reaches an acceptor in A
  # across the face, but not one in B.
  for acc_alt, across in (("A", True), ("B", False)):
    pdb_str = _LATTICE_CONTACT_PDB.replace(" O   HOH", " O  AHOH").replace(
      " O   ACA", f" O  {acc_alt}ACA")
    hier, cs = _hierarchy_and_symmetry(pdb_str)
    wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs)
    o, hs = _water_atoms(hier)
    aimed = max(_unit(h, o)[0] for h in hs.values())
    assert (aimed > 0.99) == across, (acc_alt, aimed)


def exercise_crystal_symmetry_leaves_model_fixed():
  """Symmetry equivalents join the environment as copies, so placing with a
  crystal symmetry moves no atom of the model itself."""
  hier, cs = _hierarchy_and_symmetry(_LATTICE_CONTACT_PDB)
  before = [(a.id_str(), a.xyz) for a in hier.atoms()]
  wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs)
  after = {a.id_str(): a.xyz for a in hier.atoms()}
  for atom_id, xyz in before:
    assert after[atom_id] == xyz, (
      f"{atom_id} moved from {xyz} to {after[atom_id]}")


def _min_sym_equiv_proton_gap(hier, cell_a):
  """Closest approach between a placed water proton and the lattice copy of
  any proton one cell away along a.

  The translations are applied here rather than taken from the symmetry
  machinery, so the test measures the placement instead of agreeing with it.

  Parameters
  ----------
  hier : iotbx.pdb.hierarchy.root
      Hierarchy carrying placed water H.
  cell_a : float
      Cell edge along x, in A.

  Returns
  -------
  float
      The smallest distance from a proton to a translated proton.
  """
  pts = [matrix.col(a.xyz) for ag in hier.atom_groups()
         if wp._is_water(ag.resname) for a in ag.atoms()
         if a.element.strip().upper() in ("H", "D")]
  gaps = [(p - (q + matrix.col((shift, 0.0, 0.0)))).length()
          for shift in (-cell_a, cell_a) for p in pts for q in pts]
  return min(gaps)


def exercise_sym_equiv_protons():
  """Waters donating across a lattice contact clear each other's
  symmetry-equivalent protons, not just the equivalent oxygens they aim at."""
  hier, cs = _hierarchy_and_symmetry(_CONTACT_PAIR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs)
  gap = _min_sym_equiv_proton_gap(hier, 10.0)
  assert gap >= wp._WATER_MIN_H_CLEARANCE - 1e-9, (
    f"proton sits {gap:.3f} A from a translated proton, "
    f"under the {wp._WATER_MIN_H_CLEARANCE} A clearance")


def exercise_sym_equiv_contacts_counted():
  """The clash counts include contacts with symmetry-equivalent water H.

  Places on ``_SYM_CONTACT_PDB``, whose waters already carry H, and reads
  the counts the placer reports for its state, the counts from the
  hierarchy and the program's contact listing. Each must hold the two
  contacts across symmetry once, at 1.486 and 1.600 A, the second against
  the water's own equivalent; without a crystal symmetry there are none. On
  ``_SPLIT_SYM_CONTACT_PDB`` a water split between altlocs, with an O per
  conformer, makes two contacts, each counted once.
  """
  hier, cs = _hierarchy_and_symmetry(_SYM_CONTACT_PDB)
  states = []
  wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs,
                           on_state=lambda label, stats: states.append(stats))
  assert [s[1] for s in states] == [2], states
  assert wp._water_clash_stats(hier, crystal_symmetry=cs) == states[0]
  assert wp._water_clash_stats(hier)[1] == 0
  listed = wp._worst_water_clashes(hier, crystal_symmetry=cs)
  assert approx_equal([c[0] for c in listed], [1.486, 1.600], eps=1e-6)
  assert [c[1:] for c in listed] == [
    ("HOH W 2 H1", "HOH W 3 H1 (x-1,y,z)"),
    ("HOH W 1 H1", "HOH W 1 H1 (-x,-y+1,-z+1)")], listed

  # A water split between altlocs has an O per conformer, yet each contact
  # with it counts once.
  hier, cs = _hierarchy_and_symmetry(_SPLIT_SYM_CONTACT_PDB)
  states = []
  wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs,
                           on_state=lambda label, stats: states.append(stats))
  assert [s[1] for s in states] == [2], states
  assert wp._water_clash_stats(hier, crystal_symmetry=cs) == states[0]
  listed = wp._worst_water_clashes(hier, crystal_symmetry=cs)
  assert approx_equal([c[0] for c in listed], [1.486, 1.603], eps=1e-3)
  assert [c[1:] for c in listed] == [
    ("HOH W 2 H1", "HOH W 3 H1 (A) (x-1,y,z)"),
    ("HOH W 2 H1", "HOH W 3 H1 (B) (x-1,y,z)")], listed


def exercise_sym_equiv_on_symmetry_element():
  """An atom on a symmetry element is not copied onto itself.

  Builds the symmetry environment of ``_ON_TWO_FOLD_PDB``. The two-fold
  brings the H's equivalent next to the O, so the water's residue joins the
  environment; the O, which the two-fold maps onto itself, must stay out.
  The environment holds the H's equivalent alone, traced to the H.
  """
  hier, cs = _hierarchy_and_symmetry(_ON_TWO_FOLD_PDB)
  hier.atoms().reset_i_seq()
  hiers, atoms, xyz, source = wp._symmetry_environment(
    hier, hier.atoms().extract_xyz(), cs, wp._WATER_ACCEPTOR_RADIUS)
  assert [a.name.strip() for a in atoms] == ["H1"], (
    [a.name for a in atoms], list(xyz))
  assert approx_equal(xyz[0], (-0.757, 5.586, 0.0))
  assert list(source) == [1], list(source)


def exercise_water_conformers():
  """A water split between altlocs is placed per conformer, its blank atoms
  plus one altloc's.

  Completes ``_ALTLOC_D_WATER_PDB`` stripped of its D2s: each altloc gains a
  D2 at its own occupancy. Completes ``_ALTLOC_O_WATER_PDB``: each altloc
  gains a D2 beside the shared blank D1, the water is reported once as
  carrying a single H, and its own D count as no clash. Places
  ``_WATER_SO4_ALTLOC_PDB``: the water's altloc gains two H, the sulfate's
  none. Reorients ``_ALTLOC_D_WATER_PDB``: each altloc's pair is re-placed
  at its occupancy.
  """
  def hd_by_altloc(hier):
    return {ag.altloc: sorted((a.name.strip(), a.element.strip(), a.occ)
                              for a in ag.atoms()
                              if a.element.strip() in ("H", "D"))
            for ag in hier.atom_groups()}
  pairs = {"": [], "A": [("D1", "D", 0.6), ("D2", "D", 0.6)],
           "B": [("D1", "D", 0.4), ("D2", "D", 0.4)]}

  hier = _hierarchy("\n".join(l for l in _ALTLOC_D_WATER_PDB.split("\n")
                               if " D2 " not in l))
  wp.place_water_hydrogens(hier, n_refine=0, existing_h="complete")
  assert hd_by_altloc(hier) == pairs, hd_by_altloc(hier)

  hier = _hierarchy(_ALTLOC_O_WATER_PDB)
  states = []
  res = wp.place_water_hydrogens(
    hier, n_refine=0, existing_h="complete",
    on_state=lambda label, stats: states.append(stats[1]))
  got = hd_by_altloc(hier)
  assert got == {"": [("D1", "D", 1.0)], "A": [("D2", "D", 0.5)],
                 "B": [("D2", "D", 0.5)]}, got
  assert [(p[0], p[2]) for p in res.partial_waters] == [
    ("HOH W 1", "completed")], res.partial_waters
  assert states == [0], states

  hier = _hierarchy(_WATER_SO4_ALTLOC_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  got = {ag.resname: ag.atoms_size() for ag in hier.atom_groups()}
  assert got == {"HOH": 3, "SO4": 2}, got

  hier = _hierarchy(_ALTLOC_D_WATER_PDB)
  wp.place_water_hydrogens(hier, n_refine=0, existing_h="reorient")
  assert hd_by_altloc(hier) == pairs, hd_by_altloc(hier)


def exercise_altloc_environment():
  """A water conformer ignores atoms of another altloc, its own other
  conformer included.

  Places on ``_SPLIT_WATER_PDB``, where conformer B's O blocks, and would
  itself attract, the H conformer A aims at the acceptor. Conformer A must
  still aim an H at the acceptor, the clash counts must not count A's H
  against B's, and the contact listing names an atom's altloc.
  """
  hier = _hierarchy(_SPLIT_WATER_PDB)
  states = []
  wp.place_water_hydrogens(hier, n_refine=0,
                           on_state=lambda label, stats: states.append(stats))
  ag_a = next(ag for ag in hier.atom_groups()
              if ag.resname == "HOH" and ag.altloc == "A")
  o = next(a for a in ag_a.atoms() if a.element.strip() == "O")
  hs = [a for a in ag_a.atoms() if a.element.strip() == "H"]
  acc = matrix.col((7.8, 5.0, 5.0))
  aimed = max(_unit(h, o).dot((acc - matrix.col(o.xyz)).normalize())
              for h in hs)
  assert aimed > 0.99, f"conformer A should aim at the acceptor ({aimed:.3f})"
  assert [s[1] for s in states] == [0], states
  assert wp._water_clash_stats(hier)[1] == 0
  assert wp._atom_id(hs[0]) == f"HOH W 1 {hs[0].name.strip()} (A)"


def exercise_improper_altloc_rejected():
  """An improper altloc is refused rather than read as an altloc of its own.

  Places on ``_IMPROPER_ALTLOC_PDB``, whose H1 is both blank and in altloc
  A; the placer must raise cctbx's Sorry for it.
  """
  try:
    wp.place_water_hydrogens(_hierarchy(_IMPROPER_ALTLOC_PDB), n_refine=0)
  except Sorry as e:
    assert "improper altloc" in str(e), str(e)
  else:
    raise Exception_expected


def exercise_electron_microscopy_isolated():
  """The program ignores the cell of an electron microscopy model.

  Runs the program on the lattice-contact fixture written as mmCIF, once
  recorded as X-ray diffraction and once as electron microscopy. The X-ray
  run must aim an H across the cell face at the translated acceptor; the
  electron microscopy run, whose cell is a map box rather than a lattice,
  must treat the model as isolated and aim none there.
  """
  hier, cs = _hierarchy_and_symmetry(_LATTICE_CONTACT_PDB)
  cif = hier.as_mmcif_string(crystal_symmetry=cs)
  for tag, method, across in (("xray", "X-RAY DIFFRACTION", True),
                              ("em", "ELECTRON MICROSCOPY", False)):
    text = cif.replace("\n_cell.", f"\n_exptl.method '{method}'\n_cell.", 1)
    assert "_exptl.method" in text, "fixture lacks the experiment method"
    file_name = f"tst_water_protonation_{tag}.cif"
    with open(file_name, "w") as f:
      f.write(text)
    result = run_program(program_class=wp_program.Program, logger=null_out(),
                         args=[file_name, "output.overwrite=True"])
    o, hs = _water_atoms(result.model.get_hierarchy())
    aimed = max(_unit(h, o)[0] for h in hs.values())
    if across:
      assert aimed > 0.99, f"{method}: no H aims across the face ({aimed:.3f})"
    else:
      assert aimed < 0.9, f"{method}: an H aims across the face ({aimed:.3f})"


def exercise_reorient_keeps_isotope():
  """``existing_h="reorient"`` with the automatic element re-places a water's
  protons as the element it carried, not the one its residue name implies.

  Reorients D on an HOH and H on a DOD, checking the element, the names and
  that the protons turned toward the acceptor, and an O split between
  altlocs that shares one blank D, whose altlocs each get a D pair. Then
  runs the program on the HOH carrying D at the X-ray O-H length, which must
  warn about placing D.
  """
  acc = matrix.col((2.700, 0.000, 0.000))
  for pdb_str, want in ((_BAD_DEUTERATED_PDB, "D"),
                        (_BAD_PROTONATED_PDB.replace("HOH", "DOD"), "H")):
    hier = _hierarchy(pdb_str)
    wp.place_water_hydrogens(hier, n_refine=0, existing_h="reorient")
    water = [a for a in hier.atoms() if wp._is_water(a.parent().resname)]
    o = next(a for a in water if a.element.strip() == "O")
    hs = [a for a in water if a.element.strip() in ("H", "D")]
    names = sorted(a.name.strip() for a in hs)
    assert names == [want + "1", want + "2"], names
    for a in hs:
      assert a.element.strip() == want, (a.name, a.element)
    best = max(_unit(h, o).dot((acc - matrix.col(o.xyz)).normalize())
               for h in hs)
    assert best > 0.9, f"{want} not re-placed toward the acceptor ({best:.3f})"

  hier = _hierarchy(_ALTLOC_O_WATER_PDB)
  wp.place_water_hydrogens(hier, n_refine=0, existing_h="reorient")
  for ag in hier.atom_groups():
    hs = sorted((a.name.strip(), a.element.strip()) for a in ag.atoms()
                if a.element.strip() in ("H", "D"))
    assert hs == ([] if not ag.altloc else [("D1", "D"), ("D2", "D")]), (
      ag.altloc, hs)

  file_name = "tst_water_protonation_reorient_d.pdb"
  with open(file_name, "w") as f:
    f.write(_BAD_DEUTERATED_PDB)
  log = StringIO()
  run_program(program_class=wp_program.Program, logger=log,
              args=[file_name, "existing_h=reorient", "water_geometry=xray",
                    "output.overwrite=True"])
  assert "warning: placing D" in log.getvalue(), log.getvalue()


def exercise_missing_elements_rejected():
  """A model without element columns is refused, not placed on blindly.

  Runs the placer and the program on a fixture with its element columns
  stripped. The placer must raise an AssertionError and the program a Sorry,
  both from the missing-element check, instead of skipping every water. The
  program's checks must accept an unidentified atom written as element X.
  """
  try:
    wp.place_water_hydrogens(_hierarchy(_NO_ELEMENT_PDB), n_refine=0)
  except AssertionError as e:
    assert "Uninterpretable elements" in str(e), str(e)
  else:
    raise Exception_expected

  file_name = "tst_water_protonation_no_element.pdb"
  with open(file_name, "w") as f:
    f.write(_NO_ELEMENT_PDB)
  try:
    run_program(program_class=wp_program.Program, logger=null_out(),
                args=[file_name, "output.overwrite=True"])
  except Sorry as e:
    assert "Uninterpretable elements" in str(e), str(e)
  else:
    raise Exception_expected

  dm = DataManager(["model", "phil"])
  dm.process_model_str("unknown_atom.pdb", _UNKNOWN_ATOM_PDB)
  params = CCTBXParser(program_class=wp_program.Program,
                       logger=null_out()).master_phil.extract()
  wp_program.Program(dm, params, logger=null_out()).validate()


def exercise_multi_model_rejected():
  """The program and the placer refuse a multi-model input.

  Runs the program on a two-model file and expects a Sorry naming the
  problem, raised before any placement; the placer, called directly on the
  same hierarchy, must raise an AssertionError instead of placing.
  """
  file_name = "tst_water_protonation_multi_model.pdb"
  with open(file_name, "w") as f:
    f.write(_MULTI_MODEL_PDB)
  try:
    run_program(program_class=wp_program.Program, logger=null_out(),
                args=[file_name, "output.overwrite=True"])
  except Sorry as e:
    assert "Multi-model" in str(e), str(e)
  else:
    raise Exception_expected

  try:
    wp.place_water_hydrogens(_hierarchy(_MULTI_MODEL_PDB), n_refine=0)
  except AssertionError as e:
    assert "one model" in str(e), str(e)
  else:
    raise Exception_expected


def exercise_output_file_name():
  """The output file name follows the output.* parameters.

  Runs the program with no output parameter, then with each of
  ``output.prefix``, ``output.suffix`` and ``output.serial``, and checks the
  name written (extension included).
  """
  file_name = "tst_water_protonation_naming.pdb"
  with open(file_name, "w") as f:
    f.write(_TWO_ACCEPTOR_PDB)
  for extra, want in (
      ([], "tst_water_protonation_naming_waters_protonated.pdb"),
      (["output.prefix=named"], "named_waters_protonated.pdb"),
      (["output.suffix=_h"], "tst_water_protonation_naming_h.pdb"),
      (["output.serial=2"],
       "tst_water_protonation_naming_waters_protonated_002.pdb")):
    result = run_program(program_class=wp_program.Program, logger=null_out(),
                         args=[file_name, "output.overwrite=True"] + extra)
    got = os.path.basename(result.output_file_name)
    assert got == want, (extra, got)


def run():
  """Run every exercise and print the CPU times and ``OK`` on success."""
  exercise_acceptor_directed()
  exercise_cation_repulsion()
  exercise_protonated_n_not_acceptor()
  exercise_lone_pair_directed()
  exercise_idempotent()
  exercise_reorient_existing()
  exercise_single_h_water()
  exercise_single_h_metal_annotation()
  exercise_hetatm_flag()
  exercise_heavy_cation_repulsion()
  exercise_h2_reachable_acceptor()
  exercise_h2_after_acceptor_fallback()
  exercise_refinement_reduces_clashes()
  exercise_completed_water_survives_refinement()
  exercise_element_override()
  exercise_oh_length_auto()
  exercise_environment_hydrogen_count()
  exercise_detect_neutron()
  exercise_crystal_symmetry()
  exercise_crystal_symmetry_leaves_model_fixed()
  exercise_sym_equiv_protons()
  exercise_sym_equiv_contacts_counted()
  exercise_sym_equiv_on_symmetry_element()
  exercise_water_conformers()
  exercise_altloc_environment()
  exercise_improper_altloc_rejected()
  exercise_electron_microscopy_isolated()
  exercise_reorient_keeps_isotope()
  exercise_missing_elements_rejected()
  exercise_multi_model_rejected()
  exercise_output_file_name()
  print(format_cpu_times())
  print("OK")


if __name__ == "__main__":
  run()
