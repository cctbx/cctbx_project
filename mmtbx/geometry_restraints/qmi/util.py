"""Small shared helpers for the qmi module, so cross-cutting constants and
idioms have a single home instead of drifting between modules.

Holds: the first-shell coordination cutoff, the CShM geometry-class names (the
vocabulary :mod:`shape` produces and :mod:`electronic` consumes), a chdir
context manager, and a stable residue-identity key. No heavy dependencies, so
any qmi module can import it without pulling in NumPy/xtb.
"""

from __future__ import absolute_import, division, print_function

import os
from contextlib import contextmanager


# Metal first-shell coordination cutoff (Angstrom): an atom within this of the
# metal is treated as coordinating. ~3.0 A per the protein metal-site surveys
# [9] and CheckMyMetal [10]. Shared by restraints (donor detection), qm
# (coordinating-backbone freeing), and protonation (metal-donor deprotonation)
# so the first-shell definition cannot silently diverge across the pipeline.
DEFAULT_COORDINATION_CUTOFF_A = 3.0


# common_residue_names_get_class outputs that count as a protein residue (standard
# + modified, so e.g. MSE selenomethionine qualifies). Shared by qm (coordinating
# backbone) and donors (backbone donor classification) so the two cannot diverge.
AA_RESIDUE_CLASSES = frozenset(("common_amino_acid", "modified_amino_acid"))


# CShM geometry-class names: the single vocabulary shape.py produces (as
# IDEAL_POLYHEDRA / SHAPES_BY_CN keys) and electronic.py consumes. Named here so
# a rename/typo cannot silently diverge producer from consumer (a divergence
# would mis-route the oxidation/spin prediction with no error).
LINEAR = "linear"
TRIGONAL = "trigonal"
TETRAHEDRAL = "tetrahedral"
SQUARE_PLANAR = "square_planar"
TRIGONAL_BIPYRAMIDAL = "trigonal_bipyramidal"
SQUARE_PYRAMIDAL = "square_pyramidal"
OCTAHEDRAL = "octahedral"
TRIGONAL_PRISM = "trigonal_prism"
# Higher-coordination shapes scored by shape.py but not consumed by electronic.
PENTAGONAL_BIPYRAMIDAL = "pentagonal_bipyramidal"
SQUARE_ANTIPRISM = "square_antiprism"
CUBIC = "cubic"


@contextmanager
def in_directory(path):
  """chdir into *path* for the block, restoring the prior cwd afterwards.

  xTB (via ``xtb_manager``) and ``mmtbx.reduce2`` write their scratch/output
  files into the current working directory, so a per-run scratch dir keeps
  concurrent sites from colliding and keeps the source tree clean. ``path=None``
  is a no-op (stay in the current directory).
  """
  if path is None:
    yield
    return
  saved = os.getcwd()
  os.chdir(path)
  try:
    yield
  finally:
    os.chdir(saved)


def residue_key(atom):
  """Stable residue identity for *atom*: ``(chain_id, resseq, icode, resname)``,
  all stripped so the "same residue" compares equal across modules.

  Altloc is deliberately EXCLUDED: two atoms of one residue (e.g. the two O of a
  carboxylate) must group together even when one is modelled in an alternate
  conformation, so a mixed/blank altloc must not split the residue.
  """
  ag = atom.parent()
  rg = ag.parent()
  return (rg.parent().id.strip(), rg.resseq.strip(), rg.icode.strip(),
          ag.resname.strip())
