"""Symmetry-op canonicalisation and adjacency helpers for the endo_exo
QM-region builder."""

from __future__ import absolute_import, division, print_function

from cctbx import sgtbx


def _canon_op(op):
  """Return a hash-stable copy of *op*.

  ``sgtbx.rt_mx`` has an inconsistent ``__hash__`` and ``__eq__`` after
  composition: two operations that compare equal can hash to different
  values, which breaks set and dict membership.  Round-tripping through
  the canonical xyz string makes the hash consistent.
  """
  return sgtbx.rt_mx(op.as_xyz())


def _neighbour_iseqs(adjacency, i_seq):
  """Return the set of bare neighbour i_seqs for *i_seq* in the tagged
  adjacency, dropping the per-edge ``sym_op``.

  For callers such as bond-cut detection and degree counting, to which
  only covalent connectivity matters and not the symmetry image a
  neighbour belongs to.
  """
  return {j for (j, _op) in adjacency.get(i_seq, set())}
