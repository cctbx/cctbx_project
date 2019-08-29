from __future__ import absolute_import, division, print_function
from scitbx.array_family import shared
from libtbx.test_utils import approx_equal
from six.moves import range
from six.moves import zip
try:
  from six.moves import cPickle as pickle
except ImportError:
  import pickle

def native(a):
  return [tuple(elem) for elem in a]

def exercise_stl_vector_unsigned():
  a = shared.stl_vector_unsigned()
  assert a.size() == 0
  assert len(a) == 0
  a.append([])
  assert a.size() == 1
  a.append([1])
  assert len(a) == 2
  a.append([1,2])
  assert native(a) == [(), (1,), (1, 2)]
  a.extend(a)
  a.reserve(len(a) + 100)
  assert native(a) == [(), (1,), (1, 2), (), (1,), (1, 2)]
  a.insert(1, [2,3])
  assert native(a) == [(), (2,3), (1,), (1, 2), (), (1,), (1, 2)]
  a[4] = [3,4]
  assert native(a) == [(), (2,3), (1,), (1, 2), (3,4), (1,), (1, 2)]
  a = a[2:5]
  assert native(a) == [(1,), (1, 2), (3, 4)]
  del a[0]
  assert native(a) == [(1, 2), (3, 4)]
  a.append([2,3])
  a.append([5,6])
  assert native(a) == [(1, 2), (3, 4), (2, 3), (5, 6)]
  del a[1:3]
  assert native(a) == [(1, 2), (5, 6)]
  del a[:]
  assert a.size() == 0
  a.append([1,2])
  a.append([2,3])
  a.append([3,4])
  del a[-1:]
  assert native(a) == [(1, 2), (2, 3)]
  del a[:1]
  assert native(a) == [(2, 3)]
  b = a.deep_copy()
  assert native(b) == [(2, 3)]
  a.clear()
  assert a.size() == 0
  assert b.size() == 1
  a = shared.stl_vector_unsigned([[1,2],[2,3]])
  assert native(a) == [(1, 2), (2, 3)]
  a = shared.stl_vector_unsigned(0)
  assert native(a) == []
  a = shared.stl_vector_unsigned(3)
  assert native(a) == [(),(),()]
  a = shared.stl_vector_unsigned(3, [1,2])
  assert native(a) == [(1,2),(1,2),(1,2)]

def exercise_stl_vector_double():
  a = shared.stl_vector_double()
  assert a.size() == 0
  a.append([2,3])
  assert a.size() == 1
  assert approx_equal(a[0], [2,3])
  a[0].append(4)
  assert approx_equal(a[0], [2,3,4])

def exercise_stl_set_unsigned():
  a = shared.stl_set_unsigned()
  assert a.size() == 0
  a = shared.stl_set_unsigned([(2,1)])
  assert a.size() == 1
  assert tuple(a[0]) == (1, 2)
  a = shared.stl_set_unsigned([(2,1),(3,5,2)])
  assert native(a) == [(1, 2), (2, 3, 5)]
  #
  from scitbx.array_family import flex
  b = shared.stl_vector_unsigned()
  s = flex.size_t()
  a.append_union_of_selected_arrays(arrays=b, selection=s)
  assert a.size() == 3
  assert a[2].size() == 0
  b.append([2,1])
  b.append([])
  b.append([4,3,2])
  s.append(2)
  s.append(1)
  a.append_union_of_selected_arrays(arrays=b, selection=s)
  assert tuple(a[3]) == (2,3,4)
  s[1] = 0
  a.append_union_of_selected_arrays(arrays=b, selection=s)
  assert tuple(a[4]) == (1,2,3,4)
  assert native(a) == [(1, 2), (2, 3, 5), (), (2, 3, 4), (1, 2, 3, 4)]
  #
  s = pickle.dumps(a, 1)
  l = pickle.loads(s)
  assert len(l) == len(a)
  for ae,le in zip(a,l):
    assert list(ae) == list(le)

def exercise_mat3_int():
  from scitbx.array_family import flex
  a = shared.mat3_int()
  assert a.size() == 0
  a = shared.mat3_int([list(range(9))])
  assert a.size() == 1
  assert list(a[0]) == list(range(9))

def exercise():
  exercise_stl_vector_unsigned()
  exercise_stl_vector_double()
  exercise_stl_set_unsigned()
  exercise_mat3_int()
  print("OK")

if (__name__ == "__main__"):
  exercise()
