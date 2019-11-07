from __future__ import absolute_import, division, print_function
from six.moves import zip
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
    pairs = matches.pairs()
    assert pairs.column(0).all_eq(flex.size_t_range(k.size()))
    miis_cb = pairs.column(1)
    assert perm.select(miis).all_eq(miis_cb)
  def check_ra():
    for cb_op, perm, inv_perm in zip(ra.cb_ops, ra.perms, ra.inv_perms):
      check_cb_op_perm(cb_op, perm)
      check_cb_op_perm(cb_op.inverse(), inv_perm)
  check_ra()
  assert ra.i_j_multiplication_table == [
    [0, 1, 2, 3, 4, 5],
    [1, 2, 0, 4, 5, 3],
    [2, 0, 1, 5, 3, 4],
    [3, 5, 4, 0, 2, 1],
    [4, 3, 5, 1, 0, 2],
    [5, 4, 3, 2, 1, 0]]
  assert ra.i_inv_j_multiplication_table == [
    [0, 1, 2, 3, 4, 5],
    [2, 0, 1, 5, 3, 4],
    [1, 2, 0, 4, 5, 3],
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
  from libtbx.test_utils import show_diff
  from six.moves import cStringIO as StringIO
  sio = StringIO()
  assert ra.show_summary(out=sio, prefix=": ") is ra
  assert not show_diff(sio.getvalue(), """\
: Lattice symmetry: R 3 2 :R (No. 155)
: Intensity symmetry: P 1 (No. 1)
:
: Indexing ambiguities:
:   k,l,h         3-fold    invariants:    4
:   l,h,k         3-fold    invariants:    4
:   -k,-h,-l      2-fold    invariants:    4
:   -l,-k,-h      2-fold    invariants:    4
:   -h,-l,-k      2-fold    invariants:    4
""")
  #
  ra = miller.reindexing.assistant(
    lattice_group=ms.space_group(),
    intensity_group=ms.space_group(),
    miller_indices=ra.miller_indices)
  check_ra()
  sio = StringIO()
  assert ra.show_summary(out=sio) is ra
  assert not show_diff(sio.getvalue(), """\
Lattice symmetry: R 3 2 :R (No. 155)
Intensity symmetry: R 3 2 :R (No. 155)

No indexing ambiguity.
""")
  assert ra.i_j_multiplication_table == [[0]]
  assert ra.i_inv_j_multiplication_table == [[0]]
  assert ra.i_j_inv_multiplication_table == [[0]]
  #
  ra = miller.reindexing.assistant(
    lattice_group=ms.space_group(),
    intensity_group=sgtbx.space_group_info(symbol="R 3 :R").group(),
    miller_indices=ra.miller_indices)
  check_ra()
  sio = StringIO()
  assert ra.show_summary(out=sio) is ra
  assert not show_diff(sio.getvalue(), """\
Lattice symmetry: R 3 2 :R (No. 155)
Intensity symmetry: R 3 :R (No. 146)

Indexing ambiguity:
  -h,-l,-k      2-fold    invariants:    4
""")
  assert ra.i_j_multiplication_table == [[0, 1], [1, 0]]
  assert ra.i_inv_j_multiplication_table == [[0, 1], [1, 0]]
  assert ra.i_j_inv_multiplication_table == [[0, 1], [1, 0]]
  #
  import math
  ta = math.acos(-1/3) * 180 / math.pi
  uc = uctbx.unit_cell((11,11,11,ta,ta,ta))
  ms = uc.complete_miller_set_with_lattice_symmetry(
    anomalous_flag=True,
    d_min=3)
  ra = miller.reindexing.assistant(
    lattice_group=ms.space_group(),
    intensity_group=sgtbx.space_group_info(symbol="I 4 (y+z,x+z,x+y)").group(),
    miller_indices=ms.expand_to_p1().indices())
  check_ra()
  sio = StringIO()
  assert ra.show_summary(out=sio) is ra
  assert not show_diff(sio.getvalue(), """\
Lattice symmetry: I 4 3 2 (y+z,x+z,x+y) (No. 211)
Intensity symmetry: I 4 (y+z,x+z,x+y) (No. 79)

Indexing ambiguities:
  k,l,h         3-fold    invariants:    2
  -l,-k,-h      2-fold    invariants:    4
  -h,-l,-k      2-fold    invariants:    4
  l,h,k         3-fold    invariants:    2
  -k,-h,-l      2-fold    invariants:    4
""")
  assert ra.i_j_multiplication_table == [
    [0, 1, 2, 3, 4, 5],
    [1, 4, 3, 5, 0, 2],
    [2, 5, 0, 4, 3, 1],
    [3, 2, 1, 0, 5, 4],
    [4, 0, 5, 2, 1, 3],
    [5, 3, 4, 1, 2, 0]]
  assert ra.i_inv_j_multiplication_table == [
    [0, 1, 2, 3, 4, 5],
    [4, 0, 5, 2, 1, 3],
    [2, 5, 0, 4, 3, 1],
    [3, 2, 1, 0, 5, 4],
    [1, 4, 3, 5, 0, 2],
    [5, 3, 4, 1, 2, 0]]
  assert ra.i_j_inv_multiplication_table == [
    [0, 4, 2, 3, 1, 5],
    [1, 0, 3, 5, 4, 2],
    [2, 3, 0, 4, 5, 1],
    [3, 5, 1, 0, 2, 4],
    [4, 1, 5, 2, 0, 3],
    [5, 2, 4, 1, 3, 0]]
  #
  print("OK")

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
