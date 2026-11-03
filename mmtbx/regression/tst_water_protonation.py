"""Regression tests for ``water_protonation.place_water_hydrogens``, the
H-bond-aware water-hydrogen placer, and the ``mmtbx.development.naiad``
program that wraps it.
"""

import math
import os
from io import StringIO

import iotbx.cif
import iotbx.pdb
import libtbx.load_env
from iotbx.cli_parser import CCTBXParser, run_program
from iotbx.data_manager import DataManager
from scitbx import matrix
from scitbx.array_family import flex
from libtbx.test_utils import Exception_expected, approx_equal
from libtbx.utils import Sorry, format_cpu_times, null_out

from mmtbx.hydrogens import water_protonation as wp
from mmtbx.programs import water_protonation as wp_program


# One HOH (bare O) with two acceptor O 104.5 deg apart, wider than the H-O-H
# angle: A1 at 2.6 A along +x, A2 at 2.8 A at 104.5 deg from +x (+z
# azimuth). The placer should orient H1 -> A1 (nearer) and H2 -> A2.
_TWO_ACCEPTOR_PDB = """\
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  O   ACA D   1       7.600   5.000   5.000  1.00 10.00           O
HETATM    3  O   ACB D   2       4.299   5.000   7.711  1.00 10.00           O
END
"""

# The same, without element columns, as older PDB files are written.
_NO_ELEMENT_PDB = "\n".join(l[:76].rstrip()
                            for l in _TWO_ACCEPTOR_PDB.split("\n"))

# _TWO_ACCEPTOR_PDB with its second acceptor a D atom.
_D_ELSEWHERE_PDB = _TWO_ACCEPTOR_PDB.replace(
  "O   ACB D   2       4.299   5.000   7.711  1.00 10.00           O",
  "D   ACB D   2       4.299   5.000   7.711  1.00 10.00           D")

# A water by an unidentified atom (UNX, element X) in a crystal.
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

# A water carrying a single H, with one acceptor along +x. The deposited
# proton points along +y, well off the acceptor, so completing it must build
# the partner on *that* proton's cone.
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
# off the O-H1 axis, so no point on the H-O-H cone can aim at it. ACC is
# further away and lies 104.5 deg from +x, on the cone to within a degree, so
# it is the one H2 can donate to.
_UNREACHABLE_ACCEPTOR_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   ACA D   1       2.700   0.000   0.000  1.00 10.00           O
HETATM    3  O   ACB E   1       2.425   1.400   0.000  1.00 10.00           O
HETATM    4  O   ACC F   1      -0.725  -2.807   0.000  1.00 10.00           O
END
"""

# A bare water whose H1 takes ACA (+x). On H2's cone (103.91 deg from +x) lie
# a protonated water, exactly on the cone, and ACB, 16 deg off it: H2 can
# reach both. The FAR copy moves ACB out of H2's reach.
_WATER_OR_LIGAND_ACCEPTOR_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   ACA D   1       2.600   0.000   0.000  1.00 10.00           O
HETATM    3  O   ACB E   1      -1.450   0.000   2.511  1.00 10.00           O
HETATM    4  O   HOH W   2      -0.661   2.669   0.000  1.00 10.00           O
HETATM    5  H1  HOH W   2       0.095   3.260   0.000  1.00 10.00           H
HETATM    6  H2  HOH W   2      -1.417   3.260   0.000  1.00 10.00           H
END
"""
_WATER_OR_FAR_LIGAND_ACCEPTOR_PDB = _WATER_OR_LIGAND_ACCEPTOR_PDB.replace(
  "-1.450   0.000   2.511", "-2.856   0.000   0.504")

# Six bare waters on a tight 2.6 A grid with no other acceptors: the greedy
# pass leaves close H-H contacts that refinement relaxes.
_WATER_CLUSTER_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   HOH W   2       0.000   2.600   0.000  1.00 10.00           O
HETATM    3  O   HOH W   3       2.600   0.000   0.000  1.00 10.00           O
HETATM    4  O   HOH W   4       2.600   2.600   0.000  1.00 10.00           O
HETATM    5  O   HOH W   5       5.200   0.000   0.000  1.00 10.00           O
HETATM    6  O   HOH W   6       5.200   2.600   0.000  1.00 10.00           O
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

# A bare water (W 1) with ACA 2.6 A along +x and ACB 2.9 A at 55 deg from +x,
# and a second bare water 3.4 A along -x that makes the two a cluster. H1
# takes ACA, and ACB lies 49 deg off H2's cone, beyond the donation limit.
# The even split over ACA and ACB puts each proton 24.5 deg off its acceptor.
_CLUSTER_SEARCH_PDB = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  O   ACA D   1       2.600   0.000   0.000  1.00 10.00           O
HETATM    3  O   ACB E   1       1.663   2.376   0.000  1.00 10.00           O
HETATM    4  O   HOH W   2      -3.400   0.000   0.000  1.00 10.00           O
END
"""
# The same with a third water deposited on top of the second.
_CLUSTER_SEARCH_OVERLAP_PDB = _CLUSTER_SEARCH_PDB.replace("END", (
  "HETATM    5  O   HOH W   3      -3.400   0.000   0.000  1.00 10.00"
  "           O\nEND"))


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
  """Hierarchy from an inline PDB string."""
  return iotbx.pdb.input(
    source_info=None, lines=pdb_str.split("\n")).construct_hierarchy()


def _hierarchy_and_symmetry(pdb_str):
  """``(hierarchy, crystal_symmetry)`` from an inline PDB string."""
  inp = iotbx.pdb.input(source_info=None, lines=pdb_str.split("\n"))
  return inp.construct_hierarchy(), inp.crystal_symmetry()


def _water_atoms(hier):
  """The O of the single HOH and a ``{name: atom}`` dict of its H/D."""
  water = [a for a in hier.atoms() if a.parent().resname.strip() == "HOH"]
  o = next(a for a in water if a.element.strip().upper() == "O")
  hs = {a.name.strip(): a for a in water
        if a.element.strip().upper() in ("H", "D")}
  return o, hs


def _unit(a, b):
  """Unit vector from ``b`` to ``a``, each an atom or an ``(x, y, z)``."""
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
# straight at those translates collide with each other at 1.10 A.
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

# The same water without H, and a ligand O 2.8 A away at 40 deg from the
# axis; the two-fold makes its equivalent at 40 deg on the other side, so
# the two acceptors are 80 deg apart, short of the H-O-H angle.
_TWO_FOLD_ACCEPTOR_PDB = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  O   HOH W   1       0.000   5.000   0.000  1.00 10.00           O
HETATM    2  O   EOH L   1       1.800   7.145   0.000  1.00 10.00           O
END
"""

_TWO_FOLD_BARE_PDB = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P 1 2 1
HETATM    1  O   HOH W   1       0.000   5.000   0.000  1.00 10.00           O
END
"""

# A water on an inversion centre of P-1, carrying one H along +x, and a
# ligand O 2.8 A away at 104.5 deg from it, 0.6 deg off H2's cone, at 60 deg
# azimuth from +y. Inversion maps the H onto -x, the other orientation of a
# water that cannot keep one ordered H-O-H on that site.
_ON_INVERSION_PDB = """\
CRYST1   10.000   10.000   10.000  90.00  90.00  90.00 P -1
HETATM    1  O   HOH W   1       5.000   5.000   5.000  1.00 10.00           O
HETATM    2  H1  HOH W   1       5.957   5.000   5.000  1.00 10.00           H
HETATM    3  O   EOH L   1       4.299   6.356   7.348  1.00 10.00           O
END
"""


# A water split between altlocs A and B, 1.39 A apart, and a blank acceptor
# 2.8 A from conformer A along +x. Conformer B's O sits 1.21 A from the H
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
  right O-H length and H-O-H angle. Acceptors wider apart than H-O-H leave
  both protons equally far off, in the acceptors' plane.

  Places on ``_TWO_ACCEPTOR_PDB``, then on a copy with the second acceptor
  124 deg from the first.
  """
  hier = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  atoms = list(hier.atoms())
  o, hs = _water_atoms(hier)
  assert set(hs) == {"H1", "H2"}, f"expected H1 and H2; got {sorted(hs)}"

  acc = sorted([a for a in atoms
                if a.parent().resname.strip() in ("ACA", "ACB")],
               key=lambda a: a.distance(o))
  a1, a2 = acc  # nearer first -> matches H1
  align = _unit(hs["H1"], o).dot(_unit(a1, o))
  assert align > 0.9, f"H1 should point at the nearest acceptor ({align:.3f})"
  align = _unit(hs["H2"], o).dot(_unit(a2, o))
  assert align > 0.9, (
    f"H2 should point at the second acceptor, not open space ({align:.3f})")

  non_water = [a for a in atoms if a.parent().resname.strip() != "HOH"]
  for h in (hs["H1"], hs["H2"]):
    d = min(h.distance(a) for a in non_water)
    assert d >= wp._WATER_MIN_CLEARANCE - 1e-6, (
      f"placed H too close to a non-water atom: {d:.3f} A")

  oh_target, hoh_target = wp._WATER_GEOMETRY["xray"]
  oh = (matrix.col(hs["H1"].xyz) - matrix.col(o.xyz)).length()
  assert approx_equal(oh, oh_target, eps=1e-3)
  ang = math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))
  assert approx_equal(ang, hoh_target, eps=1.0)

  # The second acceptor at 124 deg from the first: each proton 10.05 deg off.
  wide = _hierarchy(_TWO_ACCEPTOR_PDB.replace(
    "4.299   5.000   7.711", "3.434   5.000   7.321"))
  wp.place_water_hydrogens(wide, n_refine=0)
  o, hs = _water_atoms(wide)
  a1, a2 = sorted((a for a in wide.atoms()
                   if a.parent().resname.strip() in ("ACA", "ACB")),
                  key=lambda a: a.distance(o))
  off = [math.degrees(_unit(hs[n], o).angle(_unit(a, o)))
         for n, a in (("H1", a1), ("H2", a2))]
  assert max(abs(x - (124.0 - hoh_target) / 2) for x in off) < 0.05, off
  dev = max(abs(h.xyz[1] - 5.0) for h in hs.values())
  assert dev < 1e-3, f"protons should lie in the acceptors' plane ({dev:.4f})"


def exercise_cation_repulsion():
  """Both protons of a metal-coordinating water sit in the hemisphere away
  from the cation, even though the only acceptor is on the metal side.

  Places on ``_MG_WATER_PDB``.
  """
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


def exercise_protonated_water_kept():
  """A water that already has two H is left exactly as-is, as is a blank O
  carrying two D in each of altlocs A and B.

  Places on ``_PROTONATED_WATER_PDB`` and ``_ALTLOC_D_WATER_PDB`` with the
  defaults and compares every atom name and site.
  """
  for pdb_str in (_PROTONATED_WATER_PDB, _ALTLOC_D_WATER_PDB):
    hier = _hierarchy(pdb_str)
    before = [(a.name.strip(), tuple(a.xyz)) for a in hier.atoms()]
    wp.place_water_hydrogens(hier)
    after = [(a.name.strip(), tuple(a.xyz)) for a in hier.atoms()]
    assert after == before, (
      f"already-protonated water must be untouched:\n{before}\n{after}")


def exercise_protonated_n_not_acceptor():
  """An N that already carries an H is a donor, not an acceptor: H1 skips
  the nearer protonated N and points at the farther bare N instead.

  Places on ``_DONOR_N_PDB``.
  """
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
  assert best_acc > 0.9, (
    f"H should point at the bare (acceptor) N ({best_acc:.3f})")
  assert best_don < 0.9, (
    f"no H should point at the protonated (donor) N ({best_don:.3f})")


def exercise_reorient_existing():
  """``existing_h="reorient"`` strips the existing (mis-oriented) H and
  re-places them toward the acceptor; the default leaves them as-is.

  Places on ``_BAD_PROTONATED_PDB`` with the default and with
  ``existing_h="reorient"``.
  """
  acc = matrix.col((2.700, 0.000, 0.000))

  def best_align(hier):
    o, hs = _water_atoms(hier)
    return max((matrix.col(h.xyz) - matrix.col(o.xyz)).normalize().dot(
                 (acc - matrix.col(o.xyz)).normalize()) for h in hs.values())

  kept = _hierarchy(_BAD_PROTONATED_PDB)
  wp.place_water_hydrogens(kept, n_refine=0)
  align = best_align(kept)
  assert align < 0.5, (
    f"default run must not reorient existing H toward the acceptor "
    f"({align:.3f})")

  redone = _hierarchy(_BAD_PROTONATED_PDB)
  wp.place_water_hydrogens(redone, n_refine=0, existing_h="reorient")
  _, hs = _water_atoms(redone)
  assert len(hs) == 2, f"reorient must leave exactly two H; got {sorted(hs)}"
  align = best_align(redone)
  assert align > 0.9, (
    f'existing_h="reorient" should re-place an H toward the acceptor '
    f"({align:.3f})")


def exercise_single_h_water():
  """A water with one H is a possible hydroxide: untouched and reported by
  default, and completed on the *deposited* proton's cone on request.

  Places on ``_SINGLE_H_WATER_PDB`` with each ``existing_h`` mode.
  ``existing_h="reorient"`` strips and re-places it, reporting it as
  stripped.
  """
  acc = matrix.col((2.700, 0.000, 0.000))
  h1_in = matrix.col((0.000, 0.957, 0.000))

  kept = _hierarchy(_SINGLE_H_WATER_PDB)
  res = wp.place_water_hydrogens(kept, n_refine=0)
  o, hs = _water_atoms(kept)
  assert len(hs) == 1, (
    f"default must not complete a single-H water; got {sorted(hs)}")
  moved = (matrix.col(hs["H1"].xyz) - h1_in).length()
  assert moved < 1e-6, f"the deposited proton must not move ({moved:.2e} A)"
  assert [r[2] for r in res.partial_waters] == ["kept"], res.partial_waters

  done = _hierarchy(_SINGLE_H_WATER_PDB)
  res = wp.place_water_hydrogens(done, n_refine=0, existing_h="complete")
  o, hs = _water_atoms(done)
  assert len(hs) == 2, f"complete must add the partner; got {sorted(hs)}"
  moved = (matrix.col(hs["H1"].xyz) - h1_in).length()
  assert moved < 1e-6, (
    f"completing must not move the deposited proton ({moved:.2e} A)")
  O = matrix.col(o.xyz)
  oh_target, hoh_target = wp._WATER_GEOMETRY["xray"]
  ang = math.degrees(_unit(hs["H1"], o).angle(_unit(hs["H2"], o)))
  assert approx_equal(ang, hoh_target, eps=1e-3)
  oh = (matrix.col(hs["H2"].xyz) - O).length()
  assert approx_equal(oh, oh_target, eps=1e-6)
  align = _unit(hs["H2"], o).dot((acc - O).normalize())
  assert align > 0.9, (
    f"the new proton should aim at the acceptor ({align:.3f})")
  assert [r[2] for r in res.partial_waters] == ["completed"], (
    res.partial_waters)

  redone = _hierarchy(_SINGLE_H_WATER_PDB)
  res = wp.place_water_hydrogens(redone, n_refine=0, existing_h="reorient")
  o, hs = _water_atoms(redone)
  assert len(hs) == 2, f"reorient must give two H; got {sorted(hs)}"
  moved = (matrix.col(hs["H1"].xyz) - h1_in).length()
  assert moved > 0.1, (
    f"reorient must re-place the deposited proton ({moved:.3f} A)")
  assert [r[2] for r in res.partial_waters] == ["stripped"], (
    res.partial_waters)


def exercise_single_h_metal_annotation():
  """The single-H report flags a coordinating cation, using a first-shell
  cutoff rather than the looser proton-repulsion radius, and a cation that
  only crystal symmetry brings next to the water even when no water is
  placed.

  Places on ``_MG_SINGLE_H_PDB``, ``_MG_SINGLE_H_FAR_PDB`` and
  ``_MG_SINGLE_H_SYM_PDB``.
  """
  near = wp.place_water_hydrogens(
    _hierarchy(_MG_SINGLE_H_PDB), n_refine=0).partial_waters
  assert len(near) == 1, f"expected one single-H water; got {near}"
  rid, metal, action = near[0]
  assert rid == "HOH W 1", rid
  assert action == "kept", action
  assert metal is not None, (
    f"a 2.1 A Mg should be reported as coordinating; got {metal}")
  assert metal[0] == "MG", metal
  assert approx_equal(metal[1], 2.100, eps=1e-3)

  far = wp.place_water_hydrogens(
    _hierarchy(_MG_SINGLE_H_FAR_PDB), n_refine=0).partial_waters
  assert len(far) == 1 and far[0][1] is None, (
    f"a 2.9 A Mg is inside the repulsion radius but is not a bond; got {far}")

  hier, cs = _hierarchy_and_symmetry(_MG_SINGLE_H_SYM_PDB)
  sym = wp.place_water_hydrogens(
    hier, n_refine=0, crystal_symmetry=cs).partial_waters
  assert sym[0][1] is not None, sym
  assert approx_equal(sym[0][1][1], 2.100, eps=1e-3)


def exercise_hetatm_flag():
  """Placed H inherit the parent O's HETATM flag, so a HETATM water does
  not acquire ATOM-record protons.

  Places on ``_TWO_ACCEPTOR_PDB``.
  """
  hier = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  assert o.hetero, "fixture water should be a HETATM"
  for name, a in sorted(hs.items()):
    assert a.hetero == o.hetero, (
      f"{name} must match the parent O record type (hetero={a.hetero})")


def exercise_heavy_cation_repulsion():
  """The cation set is not limited to the first row: a Pt keeps both
  protons out of its hemisphere just as an Mg does.

  Places on ``_MG_WATER_PDB`` with the Mg replaced by Pt.
  """
  hier = _hierarchy(_MG_WATER_PDB.replace("MG", "PT"))
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  O = matrix.col(o.xyz)
  cd = (matrix.col((2.500, 0.000, 0.000)) - O).normalize()
  for name, h in sorted(hs.items()):
    proj = (matrix.col(h.xyz) - O).dot(cd)
    assert proj <= 1e-9, (
      f"{name} must stay out of the Pt hemisphere (proj={proj:.3f})")


def exercise_h2_reachable_acceptor():
  """H2 aims at an acceptor it can actually reach on the H-O-H cone, not at
  the nearest free one when that lies off the cone. Of the acceptors it can
  reach, H2 takes one outside the waters, and a water only when nothing
  else is in reach.

  Places on ``_UNREACHABLE_ACCEPTOR_PDB``, ``_WATER_OR_LIGAND_ACCEPTOR_PDB``
  and ``_WATER_OR_FAR_LIGAND_ACCEPTOR_PDB``.
  """
  hier = _hierarchy(_UNREACHABLE_ACCEPTOR_PDB)
  wp.place_water_hydrogens(hier, n_refine=0)
  o, hs = _water_atoms(hier)
  assert len(hs) == 2, sorted(hs)
  O = matrix.col(o.xyz)
  aca = (matrix.col((2.700, 0.000, 0.000)) - O).normalize()
  acc = (matrix.col((-0.725, -2.807, 0.000)) - O).normalize()
  d1, d2 = _unit(hs["H1"], o), _unit(hs["H2"], o)
  assert d1.dot(aca) > 0.99, (
    f"H1 should take the nearest acceptor; got {d1.dot(aca):.3f}")
  assert d2.dot(acc) > 0.9, (
    f"H2 should aim at the reachable acceptor on the cone; "
    f"alignment {d2.dot(acc):.3f}")

  for pdb_str, want in ((_WATER_OR_LIGAND_ACCEPTOR_PDB, "ACB"),
                        (_WATER_OR_FAR_LIGAND_ACCEPTOR_PDB, "HOH")):
    hier = _hierarchy(pdb_str)
    wp.place_water_hydrogens(hier, n_refine=0)
    bare = next(rg for rg in hier.residue_groups() if rg.resseq.strip() == "1"
                and rg.parent().id == "W")
    o = next(a for a in bare.atoms() if a.element.strip() == "O")
    hs = [a for a in bare.atoms() if a.element.strip() == "H"]
    assert len(hs) == 2, len(hs)
    target = next(a for a in hier.atoms() if a.parent().resname == want
                  and a.element.strip() == "O" and a.distance(o) > 0.1)
    align = max(_unit(h, o).dot(_unit(target, o)) for h in hs)
    assert align > 0.95, (
      f"a proton should aim at {want}; best alignment {align:.3f}")


def exercise_h2_after_acceptor_fallback():
  """H2 ignores H1's acceptor when the H1 fallback chose it.

  Places on ``_fallback_cage_pdb`` with the gas-phase geometry it is
  built for, where every H1 candidate clashes and the fallback settles on
  the lone acceptor. With no other acceptor, H2 must take the clash-free
  point of the sampled H-O-H cone with the most clearance.
  """
  hier = _hierarchy(_fallback_cage_pdb())
  wp.place_water_hydrogens(hier, n_refine=0, geometry="gas_phase")
  o, hs = _water_atoms(hier)
  acc = next(a for a in hier.atoms() if a.parent().resname == "ACA")
  d1 = _unit(hs["H1"], o)
  align = d1.dot(_unit(acc, o))
  assert align > 0.9999, f"H1 is not on the acceptor ({align:.6f})"
  heavy = [matrix.col(a.xyz) for a in hier.atoms()
           if a.parent().resname != "HOH"]
  def clearance(pt):
    return min((pt - x).length() for x in heavy)
  # H1 clashes, so it came from the fallback.
  h1_clearance = clearance(matrix.col(hs["H1"].xyz))
  assert h1_clearance < wp._WATER_MIN_CLEARANCE, h1_clearance
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
  a tight cluster. Basin hopping runs while a placed water clashes and stops
  once none does.

  Places on ``_WATER_CLUSTER_PDB`` without refinement, with five sweeps and
  with two basin-hopping rounds, then on ``_TWO_ACCEPTOR_PDB`` with two
  rounds. The contact listing holds the greedy pass's contacts, closest
  first.
  """
  greedy = _hierarchy(_WATER_CLUSTER_PDB)
  wp.place_water_hydrogens(greedy, n_refine=0)
  refined = _hierarchy(_WATER_CLUSTER_PDB)
  wp.place_water_hydrogens(refined, n_refine=5)

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

  # One round clears the cluster, so the second is skipped; a clash-free
  # placement runs none.
  hopped = _hierarchy(_WATER_CLUSTER_PDB)
  res = wp.place_water_hydrogens(hopped, n_refine=0, n_basin=2)
  assert res.basin_rounds_run == 1, res.basin_rounds_run
  res = wp.place_water_hydrogens(_hierarchy(_TWO_ACCEPTOR_PDB), n_refine=0,
                                 n_basin=2)
  assert res.basin_rounds_run == 0, res.basin_rounds_run


def exercise_completed_water_survives_refinement():
  """Refinement sweeps and basin kicks must hold a completed water's O-H1
  axis fixed: the deposited proton never moves and the partner stays on its
  cone.

  Completes ``_MIXED_CLUSTER_PDB`` with three sweeps and two basin-hopping
  rounds.
  """
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
  assert len(res.partial_waters) == 4, res.partial_waters
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
    assert approx_equal(ang, wp._WATER_GEOMETRY["neutron"][1], eps=0.1,
                        prefix=f"water {ag.parent().resseq.strip()}: ")


def exercise_cluster_search():
  """The cluster search turns a clustered water to the candidate whose
  protons lie closest to acceptors, here the even split over two acceptors
  that the placement, taking the nearest acceptor first, cannot reach.

  Places on ``_CLUSTER_SEARCH_PDB`` without the search, where one proton of
  W 1 points at no acceptor within the donation limit, and with it, where
  both do; then on ``_CLUSTER_SEARCH_OVERLAP_PDB``, whose two coincident
  waters give each other no direction to donate along.
  """
  for pdb_str, search, want in (
      (_CLUSTER_SEARCH_PDB, False, ["ACA"]),
      (_CLUSTER_SEARCH_PDB, True, ["ACA", "ACB"]),
      (_CLUSTER_SEARCH_OVERLAP_PDB, True, ["ACA", "ACB"])):
    hier = _hierarchy(pdb_str)
    wp.place_water_hydrogens(hier, cluster_search=search)
    bare = next(rg for rg in hier.residue_groups() if rg.resseq.strip() == "1")
    o = next(a for a in bare.atoms() if a.element.strip() == "O")
    hs = [a for a in bare.atoms() if a.element.strip() == "H"]
    assert len(hs) == 2, len(hs)
    got = []
    for a in hier.atoms():
      if a.parent().resname not in ("ACA", "ACB"):
        continue
      for h in hs:
        if _unit(h, o).dot(_unit(a, o)) >= wp._WATER_H1_ACC_ALIGN:
          got.append(a.parent().resname)
    assert sorted(got) == want, (search, sorted(got))


def exercise_cluster_search_modes():
  """Both cluster search modes keep only allowed combinations: the
  exhaustive search finds the lowest total score, iterated conditional modes
  the first local minimum from every water's current candidate.

  Two waters of two candidates each, candidate 0 the current one, whose
  second candidates may not coexist. Water 0's second candidate improves
  first and blocks water 1's, the better improvement.
  """
  scores = [flex.double([2.0, 1.0]), flex.double([2.0, 0.5])]
  allowed = {(0, 1): flex.bool([True, True, True, False])}
  assert wp._search_exhaustive(scores, allowed) == [0, 1]
  assert wp._search_icm(scores, allowed) == [1, 0]


def exercise_element_override():
  """``element="D"`` forces deuterium (named D1/D2); the default is H, or D
  when the model has D elsewhere.

  Places on ``_TWO_ACCEPTOR_PDB`` with the default and with ``element="D"``,
  then on ``_D_ELSEWHERE_PDB``.
  """
  default = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(default, n_refine=0)
  _, hs = _water_atoms(default)
  assert set(hs) == {"H1", "H2"}, f"default should place H; got {sorted(hs)}"

  deut = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(deut, n_refine=0, element="D")
  _, hs = _water_atoms(deut)
  assert set(hs) == {"D1", "D2"}, (
    f"element='D' should place D; got {sorted(hs)}")
  for d in hs.values():
    assert d.element.strip().upper() == "D", (
      f"placed atom element must be D; got {d.element}")

  # A D on another residue.
  elsewhere = _hierarchy(_D_ELSEWHERE_PDB)
  wp.place_water_hydrogens(elsewhere, n_refine=0)
  _, hs = _water_atoms(elsewhere)
  assert set(hs) == {"D1", "D2"}, f"D elsewhere should give D; got {sorted(hs)}"


def exercise_geometry_auto():
  """With no ``geometry`` the placer picks it from the model: neutron where
  D is present, X-ray otherwise, both cctbx's restraint targets for HOH,
  compared with the restraint library when it is installed. An explicit
  ``geometry`` or ``oh_length`` overrides the choice.

  Places on ``_TWO_ACCEPTOR_PDB`` and ``_D_ELSEWHERE_PDB``, with and without
  an explicit ``oh_length`` or ``geometry``, and reads the restraint
  library's HOH entry.
  """

  def oh(hier):
    o, hs = _water_atoms(hier)
    lengths = [(matrix.col(h.xyz) - matrix.col(o.xyz)).length()
               for h in hs.values()]
    assert max(lengths) - min(lengths) < 1e-6, lengths
    return lengths[0]

  def hoh(hier):
    o, hs = _water_atoms(hier)
    h1, h2 = hs.values()
    return math.degrees(_unit(h1, o).angle(_unit(h2, o)))

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
  else:
    print("chem_data missing, skipping the restraint-library check")

  # Hydrogenous model: X-ray.
  xray = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(xray, n_refine=0)
  assert approx_equal(oh(xray), xray_oh, eps=1e-6)
  assert approx_equal(hoh(xray), xray_hoh, eps=1e-6, multiplier=None)

  # A D anywhere in the model means neutron, even on another residue.
  deut = _hierarchy(_D_ELSEWHERE_PDB)
  wp.place_water_hydrogens(deut, n_refine=0)
  assert approx_equal(oh(deut), neutron_oh, eps=1e-6)

  # An explicit value wins over the heuristic.
  forced = _hierarchy(_TWO_ACCEPTOR_PDB)
  wp.place_water_hydrogens(forced, n_refine=0, oh_length=neutron_oh)
  assert approx_equal(oh(forced), neutron_oh, eps=1e-6)

  # So does an explicit geometry.
  gas = _hierarchy(_D_ELSEWHERE_PDB)
  wp.place_water_hydrogens(gas, n_refine=0, geometry="gas_phase")
  assert approx_equal((oh(gas), hoh(gas)), wp._WATER_GEOMETRY["gas_phase"],
                      eps=1e-6, multiplier=None)


def exercise_environment_hydrogen_count():
  """``count_environment_hydrogens`` counts H and D outside the waters, and
  is not fooled by a model whose only hydrogens are on its waters.

  Counts on ``_ENV_H_PDB``, on copies without the residue's H and with the
  water as DOD, and on ``_WATER_CLUSTER_PDB``.
  """
  n = wp.count_environment_hydrogens(_hierarchy(_ENV_H_PDB))
  assert n == 2, n

  # Waters keep their H; the residue loses its two.
  bare = "\n".join(l for l in _ENV_H_PDB.split("\n")
                    if " H   ASN" not in l and " D   ASN" not in l)
  hier = _hierarchy(bare)
  n = wp.count_environment_hydrogens(hier)
  assert n == 0, (
    f"hydrogens on waters must not count as environment hydrogens ({n})")
  n = sum(1 for a in hier.atoms() if a.element_is_hydrogen())
  assert n == 2, f"the water H should still be present ({n})"

  # DOD is a water alias too, so deuterium on it must not count either.
  dod = "\n".join(
    (l.replace("HOH", "DOD").replace(" H1 ", " D1 ").replace(" H2 ", " D2 ")
      .replace("           H", "           D")) if "HOH" in l else l
    for l in _ENV_H_PDB.split("\n"))
  hier = _hierarchy(dod)
  water_d = [a for a in hier.atoms()
             if a.parent().resname.strip() == "DOD"
             and a.element.strip() == "D"]
  assert len(water_d) == 2, (
    f"fixture should carry two water D; got {len(water_d)}")
  n = wp.count_environment_hydrogens(hier)
  assert n == 2, n

  # Solvent-only model.
  n = wp.count_environment_hydrogens(_hierarchy(_WATER_CLUSTER_PDB))
  assert n == 0, n


def exercise_detect_geometry():
  """``_experiment_type`` names the experiment type in plain words, and
  ``_detect_geometry`` takes the geometry it implies, electron data taking
  the neutron one, then falls back to D atoms or DOD waters.

  Reads ``_TWO_ACCEPTOR_PDB`` under three EXPDTA records, two mmCIF
  ``_exptl.method`` values and none, and two models without a record, one
  carrying a D atom and one a DOD water.
  """
  def pdb_in(lines):
    return iotbx.pdb.input(source_info=None, lines=lines.split("\n"))

  neutron = "EXPDTA    NEUTRON DIFFRACTION\n" + _TWO_ACCEPTOR_PDB
  xray = "EXPDTA    X-RAY DIFFRACTION\n" + _TWO_ACCEPTOR_PDB
  diffraction = "EXPDTA    ELECTRON DIFFRACTION\n" + _TWO_ACCEPTOR_PDB
  cif = _hierarchy(_TWO_ACCEPTOR_PDB).as_mmcif_string()
  # No experiment record, but a D atom present -> neutron by fallback.
  d_atom = """\
HETATM    1  O   HOH W   1       0.000   0.000   0.000  1.00 10.00           O
HETATM    2  D   DNR E   1       3.600   0.000   0.000  1.00 10.00           D
END
"""
  # Nor a D atom, but a DOD water -> neutron.
  dod = """\
HETATM    1  O   DOD W   1       0.000   0.000   0.000  1.00 10.00           O
END
"""

  for src, name, want in (
      (neutron, "neutron diffraction", "neutron"),
      (xray, "X-ray diffraction", "xray"),
      (diffraction, "electron diffraction", "neutron"),
      (cif + "_exptl.method 'ELECTRON MICROSCOPY'\n", "electron microscopy",
       "neutron"),
      (cif + "_exptl.method 'ELECTRON CRYSTALLOGRAPHY'\n",
       "electron diffraction", "neutron"),
      (d_atom, None, "neutron"), (dod, None, "neutron"),
      (_TWO_ACCEPTOR_PDB, None, "xray")):
    pi = pdb_in(src)
    assert wp_program._experiment_type(pi)[0] == name, (
      src.split("\n")[0], name)
    got, _ = wp_program._detect_geometry(pi, pi.construct_hierarchy())
    assert got == want, (src.split("\n")[0], got, want)


def exercise_crystal_symmetry():
  """An acceptor reachable only as a lattice translate draws an H once a
  crystal symmetry is given, and is invisible without one.

  Places on ``_LATTICE_CONTACT_PDB`` with and without its symmetry, and
  with the water and acceptor in altlocs. A translate keeps its altloc: a
  water in A reaches an acceptor in A across the face, not one in B.
  """
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
  crystal symmetry moves no atom of the model itself.

  Places on ``_LATTICE_CONTACT_PDB`` with its symmetry and compares every
  model atom's site.
  """
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

  The translations are applied here, not through the symmetry code.

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
  pts = []
  for ag in hier.atom_groups():
    if not wp._is_water(ag.resname):
      continue
    for a in ag.atoms():
      if a.element.strip().upper() in ("H", "D"):
        pts.append(matrix.col(a.xyz))
  gap = None
  for shift in (-cell_a, cell_a):
    t = matrix.col((shift, 0.0, 0.0))
    for p in pts:
      for q in pts:
        d = (p - (q + t)).length()
        if gap is None or d < gap:
          gap = d
  return gap


def exercise_sym_equiv_protons():
  """Waters donating across a lattice contact clear each other's
  symmetry-equivalent protons, not just the equivalent oxygens they aim at.

  Places on ``_CONTACT_PAIR_PDB`` with its symmetry and measures each proton
  against the lattice translates of all protons.
  """
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
  stats = wp._water_clash_stats(hier, crystal_symmetry=cs)
  assert stats == states[0], (stats, states[0])
  stats = wp._water_clash_stats(hier)
  assert stats[1] == 0, stats
  listed = wp._worst_water_clashes(hier, crystal_symmetry=cs)
  assert approx_equal([c[0] for c in listed], [1.486, 1.600], eps=1e-6)
  assert [c[1:] for c in listed] == [
    ("HOH W 2 H1", "HOH W 3 H1 [x-1,y,z]"),
    ("HOH W 1 H1", "HOH W 1 H1 [-x,-y+1,-z+1]")], listed

  # A water split between altlocs has an O per conformer, yet each contact
  # with it counts once.
  hier, cs = _hierarchy_and_symmetry(_SPLIT_SYM_CONTACT_PDB)
  states = []
  wp.place_water_hydrogens(hier, n_refine=0, crystal_symmetry=cs,
                           on_state=lambda label, stats: states.append(stats))
  assert [s[1] for s in states] == [2], states
  stats = wp._water_clash_stats(hier, crystal_symmetry=cs)
  assert stats == states[0], (stats, states[0])
  listed = wp._worst_water_clashes(hier, crystal_symmetry=cs)
  assert approx_equal([c[0] for c in listed], [1.486, 1.603], eps=1e-3)
  assert [c[1:] for c in listed] == [
    ("HOH W 2 H1", "HOH W 3 H1 altloc A [x-1,y,z]"),
    ("HOH W 2 H1", "HOH W 3 H1 altloc B [x-1,y,z]")], listed


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
    hier, cs, wp._WATER_ACCEPTOR_RADIUS)
  assert [a.name.strip() for a in atoms] == ["H1"], (
    [a.name for a in atoms], list(xyz))
  assert approx_equal(xyz[0], (-0.757, 5.586, 0.0))
  assert list(source) == [1], list(source)


def exercise_water_on_special_position():
  """A water on a special position is placed with its site symmetry.

  The two-fold of ``_ON_TWO_FOLD_PDB`` (along b through x = z = 0) maps
  (x, y, z) onto (-x, y, -z). Completing the water must put H2 along H1's
  equivalent. A bare water on that axis, with ``_TWO_FOLD_ACCEPTOR_PDB``'s
  acceptor or none, must get two H that are each other's equivalents and no
  contact with its own equivalents; with the acceptor, each H must aim
  within 13 deg of it or its equivalent, the cone passing 12 deg from both
  (half the H-O-H angle less 40 deg). On the inversion centre of
  ``_ON_INVERSION_PDB``, H1's equivalent is the water's other orientation,
  not a neighbour, so H2 must still aim at the acceptor near its cone.
  """
  hier, cs = _hierarchy_and_symmetry(_ON_TWO_FOLD_PDB)
  wp.place_water_hydrogens(hier, existing_h="complete", crystal_symmetry=cs)
  o, hs = _water_atoms(hier)
  x, y, z = hs["H1"].xyz
  off = math.degrees(_unit(hs["H2"], o).angle(_unit((-x, y, -z), o)))
  assert off < 1e-4, f"H2 is {off:.3f} deg off H1's two-fold equivalent"

  acc = (1.800, 7.145, 0.000)
  acc_equiv = (-1.800, 7.145, 0.000)
  for pdb_str in (_TWO_FOLD_ACCEPTOR_PDB, _TWO_FOLD_BARE_PDB):
    hier, cs = _hierarchy_and_symmetry(pdb_str)
    wp.place_water_hydrogens(hier, crystal_symmetry=cs)
    o, hs = _water_atoms(hier)
    x, y, z = hs["H1"].xyz
    gap = (matrix.col(hs["H2"].xyz) - matrix.col((-x, y, -z))).length()
    assert gap < 1e-6, f"H2 is {gap:.3f} A from H1's two-fold equivalent"
    stats = wp._water_clash_stats(hier, crystal_symmetry=cs)
    assert stats[1] == 0, f"contacts with own equivalents: {stats}"
    if pdb_str is _TWO_FOLD_ACCEPTOR_PDB:
      for name, h in sorted(hs.items()):
        off = min(math.degrees(_unit(h, o).angle(_unit(a, o)))
                  for a in (acc, acc_equiv))
        assert off < 13.0, f"{name} is {off:.1f} deg off the acceptors"

  hier, cs = _hierarchy_and_symmetry(_ON_INVERSION_PDB)
  wp.place_water_hydrogens(hier, existing_h="complete", crystal_symmetry=cs)
  o, hs = _water_atoms(hier)
  off = math.degrees(_unit(hs["H2"], o).angle(
    _unit((4.299, 6.356, 7.348), o)))
  assert off < 1.0, f"H2 is {off:.2f} deg off the acceptor"


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
    found = {}
    for ag in hier.atom_groups():
      hd = []
      for a in ag.atoms():
        if a.element.strip() in ("H", "D"):
          hd.append((a.name.strip(), a.element.strip(), a.occ))
      found[ag.altloc] = sorted(hd)
    return found
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
  stats = wp._water_clash_stats(hier)
  assert stats[1] == 0, stats
  atom_id = wp._atom_id(hs[0])
  assert atom_id == f"HOH W 1 {hs[0].name.strip()} altloc A", atom_id


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
  must treat the model as isolated and aim none there. Each written model
  keeps its experiment type, also when the X-ray one is written as PDB.
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
    written = iotbx.pdb.input(result.output_file_name).get_experiment_type()
    assert written.lines == [method], (method, written.lines)
  result = run_program(program_class=wp_program.Program, logger=null_out(),
                       args=["tst_water_protonation_xray.cif",
                             "output.overwrite=True",
                             "output.target_output_format=pdb"])
  assert result.output_file_name.endswith(".pdb"), result.output_file_name
  written = iotbx.pdb.input(result.output_file_name).get_experiment_type()
  assert written.lines == ["X-RAY DIFFRACTION"], written.lines


def exercise_reorient_keeps_isotope():
  """``existing_h="reorient"`` with the automatic element re-places a water's
  protons as the element it carried, not the one its residue name implies.

  Reorients D on an HOH and H on a DOD: each comes back as the same
  element, named after it and turned toward the acceptor. Reorients a bare
  HOH beside the D water, which must take D as well, and an O split between
  altlocs sharing one blank D, whose altlocs must each get a D pair. Then
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

  # A bare HOH beside the D water takes D too, though reorient strips the
  # model's only D first.
  hier = _hierarchy(_BAD_DEUTERATED_PDB.replace("END", (
    "HETATM    5  O   HOH W   2      20.000  20.000  20.000  1.00 10.00"
    "           O\nEND")))
  wp.place_water_hydrogens(hier, n_refine=0, existing_h="reorient")
  elements = sorted(a.element.strip() for a in hier.atoms()
                    if a.element.strip() in ("H", "D"))
  assert elements == ["D"] * 4, elements

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
  assert "Warning: placing D" in log.getvalue(), log.getvalue()


def exercise_missing_elements_rejected():
  """A model without element columns is refused.

  Runs the placer and the program on a fixture with its element columns
  stripped. The placer must raise an AssertionError and the program a Sorry,
  both from the missing-element check. The program's checks must accept an
  unidentified atom written as element X.
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
  same hierarchy, must raise an AssertionError.
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
  """The output file name follows the output.* parameters, and the written
  model keeps the input's experiment type.

  Runs the program on a PDB input with a joint X-ray/neutron EXPDTA record,
  with no output parameter, then with each of ``output.prefix``,
  ``output.suffix``, ``output.serial`` and
  ``output.target_output_format=mmcif``, and checks the name written
  (extension included) and the experiment type read back.
  """
  file_name = "tst_water_protonation_naming.pdb"
  with open(file_name, "w") as f:
    f.write("EXPDTA    X-RAY DIFFRACTION; NEUTRON DIFFRACTION\n"
            + _TWO_ACCEPTOR_PDB)
  for extra, want in (
      ([], "tst_water_protonation_naming_waters_protonated.pdb"),
      (["output.prefix=named"], "named_waters_protonated.pdb"),
      (["output.suffix=_h"], "tst_water_protonation_naming_h.pdb"),
      (["output.serial=2"],
       "tst_water_protonation_naming_waters_protonated_002.pdb"),
      (["output.target_output_format=mmcif"],
       "tst_water_protonation_naming_waters_protonated.cif")):
    result = run_program(program_class=wp_program.Program, logger=null_out(),
                         args=[file_name, "output.overwrite=True"] + extra)
    got = os.path.basename(result.output_file_name)
    assert got == want, (extra, got)
    written = iotbx.pdb.input(result.output_file_name).get_experiment_type()
    assert sorted(written.lines) == [
      "NEUTRON DIFFRACTION", "X-RAY DIFFRACTION"], (extra, written.lines)


def run():
  """Run every exercise and print the CPU times and ``OK`` on success."""
  exercise_acceptor_directed()
  exercise_cation_repulsion()
  exercise_protonated_n_not_acceptor()
  exercise_protonated_water_kept()
  exercise_reorient_existing()
  exercise_single_h_water()
  exercise_single_h_metal_annotation()
  exercise_hetatm_flag()
  exercise_heavy_cation_repulsion()
  exercise_h2_reachable_acceptor()
  exercise_h2_after_acceptor_fallback()
  exercise_refinement_reduces_clashes()
  exercise_completed_water_survives_refinement()
  exercise_cluster_search()
  exercise_cluster_search_modes()
  exercise_element_override()
  exercise_geometry_auto()
  exercise_environment_hydrogen_count()
  exercise_detect_geometry()
  exercise_crystal_symmetry()
  exercise_crystal_symmetry_leaves_model_fixed()
  exercise_sym_equiv_protons()
  exercise_sym_equiv_contacts_counted()
  exercise_sym_equiv_on_symmetry_element()
  exercise_water_on_special_position()
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
