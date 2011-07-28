def run(args):
  assert len(args) == 0
  from cctbx import miller
  import cctbx.miller.reindexing
  from cctbx import uctbx
  from cctbx import sgtbx
  from cctbx.array_family import flex
  uc = uctbx.unit_cell((11,11,11,81,81,81))
  ms = uc.complete_miller_set_with_lattice_symmetry(
    anomalous_flag=True,
    d_min=3)
  ra = miller.reindexing.assistant(
    lattice_group=ms.space_group(),
    intensity_group=sgtbx.space_group_info(symbol="P 1").group(),
    miller_indices=ms.expand_to_p1().indices())
  mt = flex.mersenne_twister(seed=0)
  def check_cb_op_perm(cb_op, perm):
    mi_cb = cb_op.apply(ra.miller_indices)
    miis = flex.random_permutation(size=ra.miller_indices.size())[:2]
    k = cb_op.apply(ra.miller_indices.select(miis))
    matches = miller.match_indices(k, ra.miller_indices)
    assert matches.pairs().column(0).all_eq(flex.size_t_range(k.size()))
    miis_cb = matches.pairs().column(1)
    assert perm.select(miis).all_eq(miis_cb)
  for cb_op, perm, inv_perm in zip(ra.cb_ops, ra.perms, ra.inv_perms):
    check_cb_op_perm(cb_op, perm)
    check_cb_op_perm(cb_op.inverse(), inv_perm)
  assert ra.i_j_multiplication_table == [
    [0, 1, 2, 3, 4, 5],
    [1, 2, 0, 4, 5, 3],
    [2, 0, 1, 5, 3, 4],
    [3, 5, 4, 0, 2, 1],
    [4, 3, 5, 1, 0, 2],
    [5, 4, 3, 2, 1, 0]]
  assert ra.i_j_inv_multiplication_table == [
    [0, 2, 1, 3, 4, 5],
    [1, 0, 2, 4, 5, 3],
    [2, 1, 0, 5, 3, 4],
    [3, 4, 5, 0, 2, 1],
    [4, 5, 3, 1, 0, 2],
    [5, 3, 4, 2, 1, 0]]
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
