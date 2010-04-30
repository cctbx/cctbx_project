def exercise_oset():
  from libtbx.ordered import OrderedSet as oset
  o = oset()
  assert repr(o) == "OrderedSet()"
  assert len(o) == 0
  o = oset([3,5,2,5,4,2,1])
  assert list(o) == [3, 5, 2, 4, 1]
  assert 3 in o
  assert 6 not in o
  o.add(3)
  assert len(o) == 5
  o.add(6)
  assert 6 in o
  assert list(reversed(o)) == [6,1,4,2,5,3]
  assert o.pop() == 6
  assert len(o) == 5
  assert o.pop(last=False) == 3
  assert len(o) == 4
  assert repr(o) == "OrderedSet([5, 2, 4, 1])"
  assert o == oset([5, 2, 4, 1])
  assert o != oset([5, 4, 2, 1])
  assert o == set([5, 2, 4, 1])
  assert o == set([5, 4, 2, 1])

def run(args):
  assert len(args) == 0
  exercise_oset()
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
