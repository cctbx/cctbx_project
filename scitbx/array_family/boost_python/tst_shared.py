from scitbx.array_family import shared

def exercise_std_vector_size_t():
  a = shared.std_vector_size_t()
  assert a.size() == 0
  assert len(a) == 0
  a.append([])
  assert a.size() == 1
  a.append([1])
  assert len(a) == 2
  a.append([1,2])
  assert list(a) == [(), (1,), (1, 2)]
  a.extend(a)
  assert list(a) == [(), (1,), (1, 2), (), (1,), (1, 2)]
  a.insert(1, [2,3])
  assert list(a) == [(), (2,3), (1,), (1, 2), (), (1,), (1, 2)]
  a[4] = [3,4]
  assert list(a) == [(), (2,3), (1,), (1, 2), (3,4), (1,), (1, 2)]
  a = a[2:5]
  assert list(a) == [(1,), (1, 2), (3, 4)]
  del a[0]
  assert list(a) == [(1, 2), (3, 4)]
  a.append([2,3])
  a.append([5,6])
  assert list(a) == [(1, 2), (3, 4), (2, 3), (5, 6)]
  del a[1:3]
  assert list(a) == [(1, 2), (5, 6)]
  del a[:]
  assert a.size() == 0
  a.append([1,2])
  a.append([2,3])
  a.append([3,4])
  del a[-1:]
  assert list(a) == [(1, 2), (2, 3)]
  del a[:1]
  assert list(a) == [(2, 3)]
  a.clear()
  assert a.size() == 0
  a = shared.std_vector_size_t([[1,2],[2,3]])
  assert list(a) == [(1, 2), (2, 3)]
  a = shared.std_vector_size_t(0)
  assert list(a) == []
  a = shared.std_vector_size_t(3)
  assert list(a) == [(),(),()]
  a = shared.std_vector_size_t(3, [1,2])
  assert list(a) == [(1,2),(1,2),(1,2)]

def exercise_std_set_size_t():
  a = shared.std_set_size_t()
  assert a.size() == 0
  a = shared.std_set_size_t([(2,1)])
  assert a.size() == 1
  assert a[0] == (1, 2)
  a = shared.std_set_size_t([(2,1),(3,5,2)])
  assert list(a) == [(1, 2), (2, 3, 5)]

def exercise():
  exercise_std_vector_size_t()
  exercise_std_set_size_t()
  print "OK"

if (__name__ == "__main__"):
  exercise()
