"""Adjacency helper for the endo_exo QM-region builder."""


def neighbour_iseqs(adjacency, i_seq):
  """Return the neighbours of *i_seq* in the tagged adjacency, without their
  symmetry operations.

  For callers such as bond-cut detection and degree counting, to which only
  covalent connectivity matters and not the symmetry image a neighbour
  belongs to.

  Parameters
  ----------
  adjacency : dict
      ``{i_seq: {(j_seq, rt_mx), ...}}`` as built by
      :meth:`AtomGraphBuilder.build_adjacency`.
  i_seq : int

  Returns
  -------
  set of int
  """
  return {j for (j, _op) in adjacency.get(i_seq, set())}
