from scitbx.stl import set
import pickle

def exercise_unsigned():
  s = set.unsigned()
  assert s.size() == 0
  assert len(s) == 0
  assert not 1 in s
  s.insert(1)
  assert s.size() == 1
  assert len(s) == 1
  assert 1 in s
  s.insert(1)
  assert s.size() == 1
  s.insert(3)
  assert s.size() == 2
  assert s.erase(2) == 0
  assert s.size() == 2
  assert s.erase(1) == 1
  assert s.size() == 1
  s.clear()
  assert s.size() == 0
  s = set.unsigned([4,5,8,7,5,3])
  assert s.size() == 5
  assert list(s) == [3,4,5,7,8]
  s.insert(set.unsigned(s))
  assert s.size() == 5
  s.insert([0,4,6,11])
  assert list(s) == [0,3,4,5,6,7,8,11]
  d = pickle.dumps(s)
  l = pickle.loads(d)

def exercise():
  exercise_unsigned()
  print "OK"

if (__name__ == "__main__"):
  exercise()
