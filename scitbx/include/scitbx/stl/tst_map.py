from scitbx.stl import map
import pickle

def exercise_stl_string_double():
  m = map.stl_string_double()
  assert m.size() == 0
  assert len(m) == 0
  assert not "a" in m
  assert not m.has_key("a")
  assert m.get("a", -1) == -1
  assert m.size() == 0
  assert m.setdefault("a", -2) == -2
  assert m.size() == 1
  assert len(m) == 1
  assert "a" in m
  assert m.has_key("a")
  assert m["a"] == -2
  assert m.setdefault("a", -3) == -2
  assert m.size() == 1
  assert m["a"] == -2
  m["a"] = 10
  assert m["a"] == 10
  assert m.size() == 1
  m.setdefault("b")
  assert m.size() == 2
  assert m["b"] == 0
  m["b"] = 20
  assert m.size() == 2
  assert m["b"] == 20
  del m["b"]
  assert m.size() == 1
  m["b"] = 22
  assert m.size() == 2
  assert m.erase("c") == 0
  assert m.size() == 2
  assert m.erase("a") == 1
  assert m.size() == 1
  assert m.erase("a") == 0
  assert m.size() == 1
  m.clear()
  assert m.size() == 0
  m = map.stl_string_double({})
  assert m.size() == 0
  assert m.keys() == []
  assert m.values() == []
  assert m.items() == []
  m = map.stl_string_double({"b": 2, "a": 1, "c": 3})
  assert m.size() == 3
  assert m.keys() == ["a", "b", "c"]
  assert [k for k in m] == ["a", "b", "c"]
  assert m.values() == [1, 2, 3]
  assert m.items() == zip(["a", "b", "c"], [1, 2, 3])
  d = dict(m.items())
  assert len(d) == 3
  ld = list(d)
  ld.sort()
  assert ld == list(m)
  m.update(map.stl_string_double({"x": -1, "y": -2}))
  assert m.items() == zip(["a", "b", "c", "x", "y"], [1, 2, 3, -1, -2])
  m.update({"r": 9, "s": 8})
  assert m.items() == zip(["a","b","c","r","s","x","y"], [1,2,3,9,8,-1,-2])
  assert m.popitem() == ("a", 1)
  assert m.popitem() == ("b", 2)
  d = pickle.dumps(m)
  l = pickle.loads(d)
  assert l.items() == zip(["c","r","s","x","y"], [3,9,8,-1,-2])

def exercise_stl_stl_string_stl_map_stl_string_double():
  mm = map.stl_string_stl_map_stl_string_double()
  m = mm.setdefault("a")
  assert mm["a"].size() == 0
  m["b"] = 10
  assert mm["a"].size() == 1
  m["c"] = 20
  assert mm["a"].size() == 2
  assert mm["a"]["b"] == 10
  assert mm["a"]["c"] == 20
  del mm["a"]["b"]
  assert mm["a"].size() == 1
  assert m.size() == 1

def exercise():
  exercise_stl_string_double()
  exercise_stl_stl_string_stl_map_stl_string_double()
  print "OK"

if (__name__ == "__main__"):
  exercise()
