from scitbx.array_family import shared

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

def exercise_stl_set_unsigned():
  a = shared.stl_set_unsigned()
  assert a.size() == 0
  a = shared.stl_set_unsigned([(2,1)])
  assert a.size() == 1
  assert tuple(a[0]) == (1, 2)
  a = shared.stl_set_unsigned([(2,1),(3,5,2)])
  assert [tuple(s) for s in a] == [(1, 2), (2, 3, 5)]

def exercise():
  exercise_stl_vector_unsigned()
  exercise_stl_set_unsigned()
  print "OK"

if (__name__ == "__main__"):
  exercise()
