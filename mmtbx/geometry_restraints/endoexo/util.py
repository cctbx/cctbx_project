"""Symmetry-op canonicalisation and adjacency helpers for the endoexo
QM-region builder."""

from __future__ import absolute_import, division, print_function

from cctbx import sgtbx


def _canon_op(op):
  """Return a hash-stable copy of *op*.

  cctbx's ``sgtbx.rt_mx`` has an inconsistent ``__hash__`` vs
  ``__eq__`` after composition: two operations that compare equal
  can hash to different values, which breaks set / dict membership.
  Round-tripping through the canonical xyz string fixes it.
  """
  return sgtbx.rt_mx(op.as_xyz())


def _neighbour_iseqs(adjacency, i_seq):
  """Return the set of bare neighbour i_seqs for *i_seq* in the tagged
  adjacency, dropping the per-edge ``sym_op``.

  Used by bond-cut detection and degree counting where only covalent
  connectivity matters, not which symmetry image the neighbour belongs
  to.
  """
  return {j for (j, _op) in adjacency.get(i_seq, set())}
