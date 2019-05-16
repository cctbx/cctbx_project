from __future__ import absolute_import, division, print_function
from scitbx.stl import vector
from scitbx.stl import set
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle

def exercise_unsigned():
  a = vector.unsigned()
  assert a.size() == 0
  assert len(a) == 0
  a.append(13)
  assert a.size() == 1
  a.append(12)
  assert len(a) == 2
  a.append(3)
  assert list(a) == [13, 12, 3]
  a.extend(a)
  assert list(a) == [13, 12, 3, 13, 12, 3]
  a.insert(1, 23)
  assert list(a) == [13, 23, 12, 3, 13, 12, 3]
  a[4] = 15
  assert list(a) == [13, 23, 12, 3, 15, 12, 3]
  a = a[2:5]
  assert list(a) == [12, 3, 15]
  del a[0]
  assert list(a) == [3, 15]
  a.append(23)
  a.append(56)
  assert list(a) == [3, 15, 23, 56]
  del a[1:3]
  assert list(a) == [3, 56]
  del a[:]
  assert a.size() == 0
  a.append(12)
  a.append(23)
  a.append(34)
  del a[-1:]
  assert list(a) == [12, 23]
  del a[:1]
  assert list(a) == [23]
  a.clear()
  assert a.size() == 0
  a = vector.unsigned([12,23])
  assert list(a) == [12, 23]
  a = vector.unsigned(0)
  assert list(a) == []
  a = vector.unsigned(3)
  assert list(a) == [0,0,0]
  a = vector.unsigned(3, 12)
  assert list(a) == [12,12,12]
  for elems in [(), (3,6,4)]:
    a = vector.unsigned(elems)
    d = pickle.dumps(a)
    l = pickle.loads(d)
    assert tuple(l) == elems

def exercise_set_unsigned():
  v = vector.set_unsigned()
  v.append(set.unsigned([3,7,3,4]))
  assert list(v[0]) == [3,4,7]
  v[0].insert(5)
  assert list(v[0]) == [3,4,5,7]
  v[0].clear()
  assert v[0].size() == 0

def exercise():
  exercise_unsigned()
  exercise_set_unsigned()
  print("OK")

if (__name__ == "__main__"):
  exercise()
