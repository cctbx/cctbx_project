from scitbx.array_family import flex
from scitbx.python_utils import command_line
from scitbx import matrix
from libtbx.test_utils import approx_equal, not_approx_equal
from cStringIO import StringIO
import pickle
import cPickle
import math
import time
import sys, os

def exercise_flex_grid():
  g = flex.grid()
  assert g.nd() == 0
  assert g.size_1d() == 0
  assert g.is_0_based()
  assert g.origin() == ()
  assert g.all() == ()
  assert g.last() == ()
  assert g.last(1) == ()
  assert g.last(0) == ()
  assert not g.is_padded()
  assert not g.is_trivial_1d()
  g = flex.grid((2,3,5))
  assert g.nd() == 3
  assert g.size_1d() == 30
  assert g.origin() == (0,0,0)
  assert g.all() == (2,3,5)
  assert g.last() == (2,3,5)
  assert g.last(1) == (2,3,5)
  assert g.last(0) == (1,2,4)
  assert g((0,0,0)) == 0
  assert g((1,2,4)) == 29
  assert g.is_0_based()
  assert not g.is_padded()
  assert not g.is_trivial_1d()
  assert flex.grid(1).all() == (1,)
  assert flex.grid(1,2).all() == (1,2,)
  assert flex.grid(1,2,3).all() == (1,2,3)
  assert flex.grid(1,2,3,4).all() == (1,2,3,4)
  assert flex.grid(1,2,3,4,5).all() == (1,2,3,4,5)
  assert flex.grid(1,2,3,4,5,6).all() == (1,2,3,4,5,6)
  assert flex.grid(1).set_focus(1).focus() == (1,)
  assert flex.grid(1,2).set_focus(1,2).focus() == (1,2,)
  assert flex.grid(1,2,3).set_focus(1,2,3).focus() == (1,2,3)
  assert flex.grid(1,2,3,4).set_focus(1,2,3,4).focus() == (1,2,3,4)
  assert flex.grid(1,2,3,4,5).set_focus(1,2,3,4,5).focus() == (1,2,3,4,5)
  assert flex.grid(1,2,3,4,5,6).set_focus(1,2,3,4,5,6).focus() == (1,2,3,4,5,6)
  g = flex.grid((1,2,3), (4,6,8))
  assert g.nd() == 3
  assert g.size_1d() == 60
  assert not g.is_0_based()
  assert g.origin() == (1,2,3)
  assert g.all() == (3,4,5)
  assert g.last() == (4,6,8)
  assert g.last(1) == (4,6,8)
  assert g.last(0) == (3,5,7)
  assert g((1,2,3)) == 0
  assert g((3,5,7)) == 59
  assert not g.is_valid_index((0,0,0))
  assert not g.is_padded()
  assert not g.is_trivial_1d()
  g = flex.grid((1,2,3), (4,6,8), 0)
  assert g.nd() == 3
  assert g.size_1d() == 120
  assert g.origin() == (1,2,3)
  assert g.all() == (4,5,6)
  assert g.last() == (5,7,9)
  assert g.last(1) == (5,7,9)
  assert g.last(0) == (4,6,8)
  assert g((1,2,3)) == 0
  assert g((4,6,8)) == 119
  assert not g.is_valid_index((0,0,0))
  assert not g.is_valid_index((5,0,0))
  assert not g.is_0_based()
  assert not g.is_padded()
  assert not g.is_trivial_1d()
  g.set_focus((3,-9,5))
  assert g.is_padded()
  assert g.focus() == (3,-9,5)
  assert g.focus(False) == (2,-10,4)
  g.set_focus((3,-9,5), False)
  assert g.focus() == (4,-8,6)
  assert not g.is_0_based()
  assert g.is_padded()
  assert not g.is_trivial_1d()
  s = pickle.dumps(g)
  l = pickle.loads(s)
  assert g.origin() == l.origin()
  assert g.all() == l.all()
  assert g.focus() == l.focus()
  assert g == l
  assert not g != l
  l = flex.grid((1,2,4), (4,6,8), 0).set_focus((3,-9,5))
  assert not g == l
  assert g != l
  l = flex.grid((1,2,3), (4,7,8), 0).set_focus((3,-9,5))
  assert not g == l
  assert g != l
  l = flex.grid((1,2,3), (4,6,8), 0).set_focus((4,-9,5))
  assert not g == l
  assert g != l
  g = flex.grid(1,2)
  h = flex.grid((0,0),(1,2))
  assert g == g
  assert h == h
  assert g == h
  assert h == g
  h.set_focus((1,2))
  assert not h.is_padded()
  assert g == g
  assert h == h
  assert g == h
  assert h == g
  h.set_focus((1,1))
  assert h.is_padded()
  assert h == h
  assert g != h
  assert h != g
  assert flex.grid(2).set_focus(1).is_padded()
  assert not flex.grid(1).set_focus(1).is_padded()
  assert flex.grid(1) == flex.grid(1).set_focus(1)
  assert flex.grid(1,2).set_focus(1,1).is_padded()
  assert not flex.grid(1,2).set_focus(1,2).is_padded()
  assert flex.grid(1,2) == flex.grid(1,2).set_focus(1,2)
  assert flex.grid(1,2,3).set_focus(1,2,2).is_padded()
  assert not flex.grid(1,2,3).set_focus(1,2,3).is_padded()
  assert flex.grid(1,2,3) == flex.grid(1,2,3).set_focus(1,2,3)
  assert flex.grid(1,2,3,4).set_focus(1,2,3,3).is_padded()
  assert not flex.grid(1,2,3,4).set_focus(1,2,3,4).is_padded()
  assert flex.grid(1,2,3,4) == flex.grid(1,2,3,4).set_focus(1,2,3,4)
  assert flex.grid(1,2,3,4,5).set_focus(1,2,3,4,4).is_padded()
  assert not flex.grid(1,2,3,4,5).set_focus(1,2,3,4,5).is_padded()
  assert flex.grid(1,2,3,4,5) == flex.grid(1,2,3,4,5).set_focus(1,2,3,4,5)
  g = flex.grid((1,2,3))
  assert g.shift_origin() == g
  g = flex.grid((1,2,3), (4,6,8))
  s = g.shift_origin()
  assert s.origin() == (0,0,0)
  assert s.all() == g.all()
  assert s.focus() == (3,4,5)
  assert s.focus_size_1d() == g.size_1d()
  g = flex.grid((1,2,3), (4,6,8)).set_focus((3,5,7))
  assert g.focus_size_1d() == 2*3*4
  s = g.shift_origin()
  assert s.origin() == (0,0,0)
  assert s.all() == g.all()
  assert s.focus() == (2,3,4)

def exercise_flex_constructors():
  f = flex.double()
  assert f.size() == 0
  assert f.capacity() == 0
  assert f.accessor().nd() == 1
  assert tuple(f.accessor().origin()) == (0,)
  assert tuple(f.accessor().all()) == (0,)
  assert tuple(f.accessor().last()) == (0,)
  assert tuple(f.accessor().last(1)) == (0,)
  assert tuple(f.accessor().last(0)) == (-1,)
  assert tuple(f.accessor().focus()) == (0,)
  assert tuple(f.accessor().focus(1)) == (0,)
  assert tuple(f.accessor().focus(0)) == (-1,)
  assert f.accessor().is_0_based()
  assert not f.accessor().is_padded()
  assert f.accessor().is_trivial_1d()
  assert f.nd() == 1
  assert f.is_0_based()
  assert tuple(f.origin()) == (0,)
  assert tuple(f.all()) == (0,)
  assert tuple(f.last()) == (0,)
  assert tuple(f.last(1)) == (0,)
  assert tuple(f.last(0)) == (-1,)
  assert not f.is_padded()
  assert tuple(f.focus()) == (0,)
  assert tuple(f.focus(1)) == (0,)
  assert tuple(f.focus(0)) == (-1,)
  assert f.focus_size_1d() == 0
  assert f.is_trivial_1d()
  assert tuple(f) == ()
  f = flex.double(flex.grid((2,3)))
  assert not f.is_square_matrix()
  f = flex.double(flex.grid((4,4)))
  assert f.is_square_matrix()
  f = flex.double(flex.grid((2,3,5)))
  assert f.size() == 30
  assert f.capacity() == 30
  assert list(f) == [0] * 30
  f = flex.double(flex.grid((2,3,5)), 42)
  assert f.size() == 30
  assert list(f) == [42] * 30
  f = flex.double(2)
  assert f.size() == 2
  assert tuple(f) == (0,0)
  f = flex.double(1, 42)
  assert f.size() == 1
  assert tuple(f) == (42,)
  f = flex.double((1,2,3,4,5))
  assert f.size() == 5
  assert tuple(f) == (1,2,3,4,5)
  f = flex.double([2,1,3,5])
  assert f.size() == 4
  assert tuple(f) == (2,1,3,5)
  f = flex.double(xrange(10,13))
  assert f.size() == 3
  assert tuple(f) == (10,11,12)
  assert flex.to_list(f) == [10,11,12]
  #
  for row_type in [list, tuple]:
    for column_type in [list, tuple]:
      for n_rows in range(4):
        for n_columns in range(4):
          matrix = list()
          for i_row in xrange(n_rows):
            matrix.append(column_type(flex.random_double(n_columns)))
          m = flex.double(row_type(matrix))
          if (n_rows == 0):
            assert m.focus() == (0,)
          else:
            assert m.focus() == (n_rows, n_columns)
          try: flex.double(row_type([column_type([]),column_type([1])]))
          except RuntimeError, e:
            assert str(e) == "matrix columns must have identical sizes."
          else: raise RuntimeError("Exception expected.")
          try: flex.double(row_type([column_type([0]),column_type(["x"])]))
          except TypeError, e:
            assert str(e) in [
              "bad argument type for built-in operation",
              "a float is required"]
          else: raise RuntimeError("Exception expected.")
  for arg in [[[0],""], ([0],"",)]:
    try: flex.double(arg)
    except RuntimeError, e:
      assert str(e) == \
        "argument must be a Python list or tuple of lists or tuples."
    else: raise RuntimeError("Exception expected.")

def exercise_misc():
  f = flex.double((1,2,3))
  assert f[0] == 1
  assert f[2] == 3
  assert f[-1] == 3
  assert f[-3] == 1
  assert f.front() == 1
  assert f.back() == 3
  f.fill(42)
  assert tuple(f) == (42,42,42)
  assert f.capacity() == 3
  f.reserve(10)
  assert f.capacity() == 10
  assert f.size() == 3
  f.reserve(5)
  assert f.capacity() == 10
  assert f.size() == 3
  f = flex.double((1,2,3,4,5,6))
  fc = f.deep_copy()
  assert fc.id() != f.id()
  assert tuple(fc) == (1,2,3,4,5,6)
  fc = f.shallow_copy()
  assert fc.id() == f.id()
  assert tuple(fc) == (1,2,3,4,5,6)
  fc.resize(flex.grid((2,3)))
  assert f.nd() == 1
  assert fc.nd() == 2
  f = flex.double((1,2,3))
  g = flex.double(f)
  assert g.id() == f.id()
  f1 = f.as_1d()
  assert f1.id() == f.id()
  assert tuple(f1) == (1,2,3)
  f = flex.double(flex.grid((2,3)))
  assert f.nd() == 2
  f1 = f.as_1d()
  assert f1.id() == f.id()
  assert f1.nd() == 1
  assert f1.size() == 6
  g = flex.grid((1,2,3), (4,6,8)).set_focus((3,5,7))
  f = flex.double(g)
  assert f.focus_size_1d() == 2*3*4
  assert f.shift_origin().accessor() == g.shift_origin()
  b = flex.bool([0,0,1,0,1,1,1,0,0,1,1,0,0,0])
  assert b.md5().hexdigest() == "a3a1ff7423c672e6252003c23ad5420f"
  b = flex.bool([0,0,1,0,1,0,1,0,0,1,1,0,0,0])
  assert b.md5().hexdigest() == "bc115dabbd6dc87323302b082152be14"
  #
  class old_style:
    def __init__(self, elems):
      self.elems = tuple(elems)
    def __len__(self):
      return len(self.elems)
    def __getitem__(self, i):
      return self.elems[i]
  class new_style(object):
    def __init__(self, elems):
      self.elems = tuple(elems)
    def __len__(self):
      return len(self.elems)
    def __getitem__(self, i):
      return self.elems[i]
  assert flex.double(old_style([2,3,4])).all_eq(flex.double([2,3,4]))
  assert flex.double(new_style([3,4,5])).all_eq(flex.double([3,4,5]))
  for s in ["", u""]:
    try: flex.double(s)
    except Exception, e:
      assert str(e).startswith("Python argument types in")
    else: raise RuntimeError("Exception expected.")
  #
  a = flex.double(12)
  assert a.focus() == (12,)
  a.reshape(flex.grid(3,4))
  assert a.focus() == (3,4)
  a.reshape(flex.grid(2,3,2))
  assert a.focus() == (2,3,2)
  try: a.reshape(flex.grid(5,6))
  except RuntimeError, e:
    assert str(e).find("SCITBX_ASSERT(grid.size_1d() == a.size())") > 0
  else: raise RuntimeError("Exception expected.")

def exercise_1d_slicing_core(a):
  assert tuple(a[:]) == (1,2,3,4,5)
  assert tuple(a[::]) == (1,2,3,4,5)
  assert tuple(a[0:]) == (1,2,3,4,5)
  assert tuple(a[1:]) == (2,3,4,5)
  assert tuple(a[-2:]) == (4,5)
  assert tuple(a[-1:]) == (5,)
  assert tuple(a[-6::]) == (1,2,3,4,5) # Numeric-21.0 a[-6:] is different
  assert tuple(a[-7::]) == (1,2,3,4,5) # Numeric-21.0 a[-7:] is different
  assert tuple(a[:-1]) == (1,2,3,4)
  assert tuple(a[:-2]) == (1,2,3)
  assert tuple(a[:-4]) == (1,)
  assert tuple(a[:-5]) == ()
  assert tuple(a[:-6:]) == () # Numeric-21.0 a[:-6] is different
  assert tuple(a[:-7:]) == () # Numeric-21.0 a[:-7] is different
  assert tuple(a[:0]) == ()
  assert tuple(a[:1]) == (1,)
  assert tuple(a[:2]) == (1,2)
  assert tuple(a[:3]) == (1,2,3)
  assert tuple(a[:4]) == (1,2,3,4)
  assert tuple(a[:5]) == (1,2,3,4,5)
  assert tuple(a[:6]) == (1,2,3,4,5)
  assert tuple(a[::1]) == (1,2,3,4,5)
  assert tuple(a[::2]) == (1,3,5)
  assert tuple(a[1::2]) == (2,4)
  assert tuple(a[1:-1:2]) == (2,4)
  assert tuple(a[1:-2:2]) == (2,)
  assert tuple(a[4:2]) == ()
  assert tuple(a[4:2:-1]) == (5,4)
  assert tuple(a[2:4:-1]) == ()
  assert tuple(a[::-1]) == (5,4,3,2,1)
  assert tuple(a[::-2]) == (5,3,1)
  assert tuple(a[-1::-2]) == (5,3,1)
  assert tuple(a[-1:1:-2]) == (5,3)
  assert tuple(a[-1:2:-2]) == (5,)
  assert tuple(a[3:3:0]) == ()

def exercise_1d_slicing():
  exercise_1d_slicing_core(flex.int((1,2,3,4,5)))
  try:
    import Numeric
  except ImportError:
    pass
  else:
    print "Testing compatibility with Numeric slicing...",
    exercise_1d_slicing_core(Numeric.array((1,2,3,4,5)))
    print "OK"
  assert list(flex.slice_indices(5, slice(10))) == [0,1,2,3,4]
  assert list(flex.slice_indices(5, slice(0))) == []
  assert list(flex.slice_indices(5, slice(0,10,2))) == [0,2,4]
  assert list(flex.slice_indices(5, slice(1,10,2))) == [1,3]
  assert list(flex.slice_indices(5, slice(1,3,2))) == [1]
  assert list(flex.slice_indices(5, slice(4,0,-2))) == [4,2]

def exercise_push_back_etc():
  a = flex.double(3)
  assert a.size() == 3
  assert tuple(a) == (0, 0, 0)
  a = flex.double()
  assert a.size() == 0
  a.assign(3, 1)
  assert tuple(a) == (1, 1, 1)
  a.push_back(2)
  assert tuple(a) == (1, 1, 1, 2)
  a.insert(0, 3)
  assert tuple(a) == (3, 1, 1, 1, 2)
  a.insert(2, 3, 4)
  assert tuple(a) == (3, 1, 4, 4, 4, 1, 1, 2)
  a.pop_back()
  assert tuple(a) == (3, 1, 4, 4, 4, 1, 1)
  a.erase(2)
  assert tuple(a) == (3, 1, 4, 4, 1, 1)
  a.erase(3, 5)
  assert tuple(a) == (3, 1, 4, 1)
  a.resize(6)
  assert tuple(a) == (3, 1, 4, 1, 0, 0)
  a.resize(8, -1)
  assert tuple(a) == (3, 1, 4, 1, 0, 0, -1, -1)
  a.clear()
  assert a.size() == 0
  a = flex.double((0, 1, 2, 3))
  b = flex.double((4, 5, 6))
  a.extend(b)
  assert tuple(a) == (0, 1, 2, 3, 4, 5, 6)
  g = flex.grid((2,3,5))
  f = flex.double(g, 11)
  g = flex.grid((2,3,7))
  f.resize(g)
  assert list(f) == [11] * 30 + [0] * 12
  g = flex.grid((2,3,9))
  f.resize(g, 22)
  assert list(f) == [11] * 30 + [0] * 12 + [22] * 12
  g = flex.grid((1,2,3))
  f.resize(g)
  assert list(f) == [11] * 6
  b = flex.double([10,13,17])
  assert list(a.reversed()) == [6,5,4,3,2,1,0]
  assert list(a.reversed().reversed()) == range(7)
  assert list(a.concatenate(b)) == [0,1,2,3,4,5,6,10,13,17]
  assert list(b.concatenate(a)) == [10,13,17,0,1,2,3,4,5,6]

def exercise_setitem():
  a = flex.double(2)
  a[0] = 11
  a[1] = 12
  assert tuple(a) == (11, 12)
  g = flex.grid((2,3))
  a = flex.double(g)
  for i in xrange(2):
    for j in xrange(3):
      a[(i,j)] = i * 3 + j
  assert list(a) == range(6)

def exercise_select():
  from scitbx import stl
  import scitbx.stl.vector
  a = flex.double((1,2,3,4,5))
  b = flex.bool((0,1,0,1,1))
  c = flex.size_t((0,2))
  d = stl.vector.unsigned((1,3))
  assert tuple(a.select(b)) == (2,4,5)
  assert tuple(a.select(indices=c)) == (1,3)
  assert tuple(a.select(indices=c, reverse=False)) == (1,3)
  p = flex.size_t([0,3,1,2,4])
  assert list(a.select(p)) == [1,4,2,3,5]
  assert list(a.select(p).select(p, reverse=True)) == list(a)
  for i_trial in xrange(10):
    p = flex.random_permutation(size=a.size())
    assert list(a.select(p).select(p, reverse=True)) == list(a)
  assert tuple(a.set_selected(b, flex.double((7,8,9)))) == (1,7,3,8,9)
  assert tuple(a.set_selected(c, flex.double((-1,-2)))) == (-1,7,-2,8,9)
  assert tuple(a.set_selected(d, flex.double((-1,-2)))) == (-1,-1,-2,-2,9)
  assert tuple(a.set_selected(flex.bool(5,0), flex.double()))==(-1,-1,-2,-2,9)
  assert tuple(a.set_selected(flex.size_t(), flex.double()))==(-1,-1,-2,-2,9)
  assert tuple(a.set_selected(stl.vector.unsigned(), flex.double())) \
      == (-1,-1,-2,-2,9)
  assert tuple(a.set_selected(b, -4)) == (-1,-4,-2,-4,-4)
  assert tuple(a.set_selected(c, -3)) == (-3,-4,-3,-4,-4)
  assert tuple(a.set_selected(d, -2)) == (-3,-2,-3,-2,-4)
  assert tuple(a.set_selected(flex.bool(5, 0), -9)) == (-3,-2,-3,-2,-4)
  assert tuple(a.set_selected(flex.size_t(), -9)) == (-3,-2,-3,-2,-4)
  assert tuple(a.set_selected(stl.vector.unsigned(), -9)) == (-3,-2,-3,-2,-4)
  for i,v in enumerate([1,2,3,4]):
    a = flex.double([1,2,3,4])
    a.resize(flex.grid(2,2))
    b = a.deep_copy()
    b[i] *= 10
    assert list(a.set_selected(a==v, v*10)) == list(b)
  #
  a = flex.double([1,-2,3])
  i = flex.size_t([0,2])
  v = flex.double([6,-4])
  assert a.add_selected(indices=i, values=v) is a
  assert approx_equal(a, [7,-2,-1])
  #
  a = flex.double((1,2,3,4,5))
  b = flex.size_t((3,1,0,4,2))
  assert tuple(a.select(b)) == (4,2,1,5,3)
  b = flex.size_t((1,4,2))
  assert tuple(a.select(b)) == (2,5,3)
  b = flex.size_t((2,4,1,2,4))
  assert tuple(a.select(b)) == (3,5,2,3,5)
  a = flex.size_t((1,2,3))
  for i in xrange(3):
    for expected in ([1, 2, 3],
                     [1, 3, 2],
                     [2, 1, 3],
                     [2, 3, 1],
                     [3, 1, 2],
                     [3, 2, 1]):
      assert list(a) == expected
      assert a.next_permutation() == (expected != [3, 2, 1])
  a = flex.bool((0,1,0,1,1))
  assert tuple(a.as_int()) == (0,1,0,1,1)
  assert tuple(a.as_double()) == (0,1,0,1,1)
  assert tuple(a.iselection()) == (1,3,4)
  assert tuple(a.iselection(test_value=True)) == (1,3,4)
  assert tuple(a.iselection(test_value=False)) == (0,2)
  a = flex.bool((0,0,0,0,0))
  assert tuple(a.iselection(False)) == (0,1,2,3,4)
  assert tuple(a.iselection(True)) == ()
  isel = stl.vector.unsigned([1,3])
  a = flex.bool(size=5, iselection=isel)
  assert tuple(a) == (False, True, False, True, False)
  isel2 = stl.vector.unsigned([1,4])
  assert tuple(flex.union(size=5, iselections=[isel,isel2])) \
      == (False, True, False, True, True)
  assert tuple(flex.intersection(size=5, iselections=[isel,isel2])) \
      == (False, True, False, False, False)
  assert list(flex.intersection(size=5,iselections=[isel,isel2]).iselection())\
      == [1]
  isel = flex.size_t([1,4])
  a = flex.bool(size=5, iselection=isel)
  assert tuple(a) == (False, True, False, False, True)
  isel2 = flex.size_t([0,4])
  assert tuple(flex.union(size=5, iselections=[isel,isel2])) \
      == (True, True, False, False, True)
  assert tuple(flex.intersection(size=5, iselections=[isel,isel2])) \
      == (False, False, False, False, True)
  assert list(flex.intersection(size=5,iselections=[isel,isel2]).iselection())\
      == [4]
  #
  a = flex.bool([True,False,False,True,True])
  b = flex.size_t([0,1,2,3])
  assert list(a.filter_indices(b)) == [0,3]
  #
  def iselection_intersection(a, b):
    return list(flex.size_t(a).intersection(other=flex.size_t(b)))
  assert iselection_intersection([], []) == []
  assert iselection_intersection([1,2,3,4], [1,2,3,4]) == [1,2,3,4]
  assert iselection_intersection([], [1,2,3,4]) == []
  assert iselection_intersection([1,2,3,4], []) == []
  assert iselection_intersection([1], [1,2,3,4]) == [1]
  assert iselection_intersection([1,2,3,4], [1]) == [1]
  assert iselection_intersection([4], [1,2,3,4]) == [4]
  assert iselection_intersection([1,2,3,4], [4]) == [4]
  assert iselection_intersection([1,4], [1,2,3,4]) == [1,4]
  assert iselection_intersection([1,2,3,4], [1,4]) == [1,4]
  assert iselection_intersection([1,4,5], [1,2,3,4]) == [1,4]
  assert iselection_intersection([1,2,3,4], [1,4,5]) == [1,4]
  assert iselection_intersection([1,2,3,4], [2]) == [2]
  assert iselection_intersection([2], [1,2,3,4]) == [2]
  assert iselection_intersection([1,2,3,4], [2,3]) == [2,3]
  assert iselection_intersection([2,3], [1,2,3,4]) == [2,3]
  assert iselection_intersection([1,2,3,4], [2,4]) == [2,4]
  assert iselection_intersection([2,4], [1,2,3,4]) == [2,4]

def exercise_from_stl_vector():
  from scitbx import stl
  import scitbx.stl.vector
  assert list(flex.size_t(stl.vector.unsigned([2,5,9]))) == [2,5,9]
  assert list(flex.double(stl.vector.double([3,-6,10]))) == [3,-6,10]

def exercise_operators():
  a = flex.bool((0, 1, 0, 1))
  b = flex.bool((0, 1, 1, 0))
  assert tuple(~a) == (1, 0, 1, 0)
  assert tuple(a & b) == (0, 1, 0, 0)
  assert tuple(a | b) == (0, 1, 1, 1)
  a = flex.int((4, 9))
  b = flex.int((2, 3))
  assert tuple(-a) == (-4, -9)
  assert tuple(a + b) == (6, 12)
  assert tuple(a - b) == (2, 6)
  assert tuple(a * b) == (8, 27)
  assert tuple(a / b) == (2, 3)
  assert tuple(a % b) == (0, 0)
  assert tuple(a + 3) == (7, 12)
  assert tuple(3 + a) == (7, 12)
  assert tuple(a - 4) == (0, 5)
  assert tuple(4 - a) == (0, -5)
  assert tuple(a * 5) == (20, 45)
  assert tuple(5 * a) == (20, 45)
  assert tuple(a / 2) == (2, 4)
  assert tuple(9 / a) == (2, 1)
  assert tuple(a % 2) == (0, 1)
  assert tuple(13 % a) == (1, 4)
  assert flex.sum(a) == 13
  a = flex.int((4, 9))
  b = flex.int((2, 12))
  assert tuple(a == b) == (0, 0)
  assert tuple(a != b) == (1, 1)
  assert tuple(a < b) == (0, 1)
  assert tuple(a > b) == (1, 0)
  assert tuple(a <= b) == (0, 1)
  assert tuple(a >= b) == (1, 0)
  assert tuple(a == 9) == (0, 1)
  assert tuple(a.as_double()) == (4, 9)
  assert a.all_eq(a)
  assert not a.all_eq(b)
  assert not a.all_eq(4)
  assert a.all_ne(b)
  assert not a.all_ne(a)
  assert not a.all_ne(4)
  assert a.all_ne(5)
  assert not a.all_lt(b)
  assert a.all_lt(10)
  assert not a.all_gt(b)
  assert a.all_gt(3)
  assert not a.all_le(b)
  assert a.all_le(9)
  assert not a.all_ge(b)
  assert a.all_ge(2)
  assert flex.order(a, b) == cmp(tuple(a), tuple(b))
  assert flex.order(b, a) == cmp(tuple(b), tuple(a))
  assert flex.order(a, a) == 0
  b = a.deep_copy()
  assert flex.order(a, b) == 0
  #
  a = flex.double([1,2,3])
  assert approx_equal(a*(2+3j), [(2+3j), (4+6j), (6+9j)])
  assert approx_equal((3+2j)*a, [(2j+3), (4j+6), (6j+9)])
  assert approx_equal(a*3, [3,6,9])
  assert approx_equal(4*a, [4,8,12])

def exercise_bool_inplace_operators():
  a = flex.bool((0, 1, 0, 1))
  b = flex.bool((0, 1, 1, 0))
  a &= b
  assert tuple(a) == (0, 1, 0, 0)
  a |= flex.bool((1, 0, 1, 0))
  assert tuple(a) == (1, 1, 1, 0)
  assert a.count(0) == 1
  assert a.count(1) == 3
  a &= 1
  assert tuple(a) == (1, 1, 1, 0)
  a &= 0
  assert tuple(a) == (0, 0, 0, 0)
  a |= 1
  assert tuple(a) == (1, 1, 1, 1)
  a |= 0
  assert tuple(a) == (1, 1, 1, 1)

def exercise_arith_inplace_operators():
  a = flex.int((4, 9))
  a += 3
  assert tuple(a) == (7, 12)
  a -= 3
  assert tuple(a) == (4, 9)
  a *= 3
  assert tuple(a) == (12, 27)
  a /= 3
  assert tuple(a) == (4, 9)
  a %= 3
  assert tuple(a) == (1, 0)
  a += 5
  a += a
  assert tuple(a) == (12, 10)
  a -= flex.int((4, 3))
  assert tuple(a) == (8, 7)
  if ("".join([str(n) for n in sys.version_info[:3]]) > "221"):
    a *= a
  else:
    a = flex.int((64, 49))
  assert tuple(a) == (64, 49)
  a /= flex.int((2, 1))
  assert tuple(a) == (32, 49)
  a %= flex.int((15, 14))
  assert tuple(a) == (2, 7)

def exercise_functions():
  a = flex.int((-1, 0, 1))
  assert tuple(flex.abs(a)) == (1, 0, 1)
  assert tuple(flex.pow2(a)) == (1, 0, 1)
  assert a.count(0) == 1
  assert a.count(2) == 0
  a = flex.int((1,1,1,3,3,3,3))
  assert a.count(1) == 3
  assert a.count(3) == 4
  a = flex.double((1, 0, 3, 2))
  b = flex.double((4, 5, 6))
  assert a.count(3) == 1
  assert b.count(3) == 0
  assert flex.min_index(a) == 1
  assert flex.min_index(b) == 0
  assert flex.max_index(a) == 2
  assert flex.max_index(b) == 2
  assert flex.min(a) == 0
  assert flex.min(b) == 4
  assert flex.max(a) == 3
  assert flex.max(b) == 6
  assert flex.sum(a) == 6
  assert flex.sum(b) == 15
  assert approx_equal(a.norm(), 3.74165738677)
  assert approx_equal(b.norm(), 8.77496438739)
  assert flex.product(a) == 0
  assert flex.product(b) == 120
  assert flex.mean(a) == 6. / 4
  assert flex.mean(b) == 15. / 3
  assert flex.mean_sq(a) == (1+3.*3.+2.*2.) / 4
  assert approx_equal(flex.mean_sq(b), (4.*4.+5.*5.+6.*6.) / 3)
  a = flex.double((-2, 0, 3))
  assert tuple(flex.pow(a, 2)) == (4, 0, 9)
  a = flex.double((2, 0, 3))
  assert list(flex.sqrt(a)) == [math.sqrt(x) for x in a]
  b = flex.double((1, 1, 1))
  assert (flex.mean(a) - flex.mean_weighted(a, b)) < 1.e-6
  assert (flex.mean_sq(a) - flex.mean_sq_weighted(a, b)) < 1.e-6
  assert a.all_approx_equal(a)
  assert a.all_approx_equal(other=a, tolerance=1.e-6)
  assert not a.all_approx_equal(other=b)
  assert a.all_approx_equal(b, tolerance=3)
  assert not a.all_approx_equal(1)
  assert b.all_approx_equal(other=1)
  assert a.all_approx_equal(other=1, tolerance=3)
  #
  a = flex.double([-3.67,-0.123,-0.678,0.321,0.765,8.01])
  b = a.round()
  assert approx_equal(b, [-4,0,-1,0,1,8])
  for n_digits in [-2,1,0,1,2]:
    b = a.round(n_digits=n_digits)
    for x,y in zip(a,b):
      assert approx_equal(round(x, n_digits), y)
  #
  a = flex.std_string(["a"]*3 + ["b"]*4)
  assert a.count("a") == 3
  assert a.count("b") == 4
  assert a.count("c") == 0
  list(flex.split_lines(multi_line_string="")) == []
  for multi_line_string in [
    "",
    "\n",
    "\r",
    "\r\n",
    "\n\r",
    "a\nb\n",
    "a\nb",
    "a\nb",
    "a\r\nb\n",
    "a\r\n\nb\n",
    "\ra\n\r\nb\n",
    "\ra\r\n\nb\n"]:
    for keep_ends in [False, True]:
      for count_lines_first in [True, False]:
        assert list(flex.split_lines(
                      multi_line_string,
                      keep_ends=keep_ends,
                      count_lines_first=count_lines_first)) \
            == multi_line_string.splitlines(keep_ends)
  #
  a = flex.int(range(-2,2) + range(2,4) + range(1,3))
  assert a.counts().items() == [(-2,1),(-1,1),(0,1),(1,2),(2,2),(3,1)]
  assert a.counts(max_keys=6).items() == a.counts().items()
  a = flex.long(range(-2,2) + range(-1,3) + range(1,3))
  assert a.counts().items() == [(-2,1),(-1,2),(0,2),(1,3),(2,2)]
  assert a.counts(max_keys=5).items() == a.counts().items()
  a = flex.size_t(range(1,2) + range(1,3) + range(1,4))
  assert a.counts().items() == [(1, 3), (2, 2), (3, 1)]
  assert a.counts(max_keys=3).items() == a.counts().items()
  try: a.counts(max_keys=2)
  except RuntimeError, e:
    assert str(e) == "scitbx::af::counts::limited: max_keys exceeded."
  else: raise RuntimeError("Exception expected.")
  #
  x = flex.double([-6.3,7.2])
  assert approx_equal(flex.fmod(x, 5), [-1.3, 2.2])
  assert approx_equal(flex.fmod_positive(x, 5), [3.7, 2.2])

def exercise_complex_functions():
  c = 1+2j
  x = flex.complex_double((c,))
  y = flex.real(x)
  assert approx_equal(y[0], 1)
  y = flex.imag(x)
  assert approx_equal(y[0], 2)
  y = flex.conj(x)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, -c.imag)
  a = flex.abs(x)
  assert approx_equal(a[0], abs(c))
  p = flex.arg(x)
  y = flex.polar(a, p)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  p = flex.arg(x, 0)
  y = flex.polar(a, p, 0)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  p = flex.arg(x, 1)
  y = flex.polar(a, p, 1)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  y = flex.polar(a, p, 0)
  d = y[0]
  assert not_approx_equal(d.real, c.real)
  assert not_approx_equal(d.imag, c.imag)
  p = flex.arg(x, 0)
  y = flex.polar(x, p)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  y = flex.polar(x, p, 0)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  y = flex.polar(x, p, 1)
  d = y[0]
  assert not_approx_equal(d.real, c.real)
  assert not_approx_equal(d.imag, c.imag)
  y = flex.polar(a, x)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  y = flex.polar(x, x)
  d = y[0]
  assert approx_equal(d.real, c.real)
  assert approx_equal(d.imag, c.imag)
  y = flex.polar(1, p)
  d = y[0]
  assert approx_equal(abs(d), 1)
  assert approx_equal(flex.arg(y)[0], p[0])
  y = flex.polar(a, math.pi/2)
  d = y[0]
  assert approx_equal(d.real, 0)
  assert approx_equal(d.imag, a[0])
  y = flex.double([2])
  assert approx_equal(x*y, [(2+4j)])
  assert approx_equal(y*x, [(2+4j)])

def exercise_sort():
  for flex_type in (flex.int, flex.size_t, flex.double):
    x = flex_type((3,1,2))
    p = flex.sort_permutation(data=x)
    assert tuple(p) == (1,2,0)
    assert approx_equal(x.select(p), (1,2,3))
    p = flex.sort_permutation(data=x, reverse=False)
    assert tuple(p) == (1,2,0)
    assert approx_equal(x.select(p), (1,2,3))
    p = flex.sort_permutation(x, True)
    assert tuple(p) == (0,2,1)
    assert approx_equal(x.select(p), (3,2,1))
  for i_trial in xrange(10):
    a = flex.size_t([0,0,0,1,1,2,2,2,3,4,4])
    if (i_trial):
      a = a.select(flex.sort_permutation(flex.random_double(size=a.size())))
    x = flex.random_double(size=5)
    p = flex.sort_permutation(data=x)
    pp = p.inverse_permutation()
    assert pp.inverse_permutation().all_eq(p)
    ap = pp.select(a)
    xp = x.select(p)
    for i,j in zip(a,ap):
      assert x[i] == xp[j]

def exercise_random():
  mt = flex.mersenne_twister()
  assert mt.random_size_t_min() == 0
  assert mt.random_size_t_max() == 4294967295
  assert mt.random_size_t() == 1791095845
  assert approx_equal(mt.random_double(), 0.997184808365)
  for i in xrange(3):
    for j in xrange(2):
      if (j == 0): mt = flex.mersenne_twister()
      else: mt.seed()
      assert tuple(mt.random_size_t(3)) \
          == (1791095845, 4282876139L, 3093770124L)
      assert approx_equal(mt.random_double(3),
       (0.9325573593386588, 0.12812444792935673, 0.99904051532414473))
      if (j == 0): mt = flex.mersenne_twister(seed=4357)
      else: mt.seed(value=4357)
      assert tuple(mt.random_size_t(size=3)) \
          == (2983900864L, 1547366158, 1775641839)
      assert approx_equal(mt.random_double(size=3),
        (0.10064729869939604, 0.89184217257908471, 0.20721445761797463))
  assert flex.random_size_t(size=3).size() == 3
  assert flex.random_double(size=3).size() == 3
  for i_trial in xrange(10):
    a = flex.random_size_t(size=100000, modulus=10)
    assert a.size() == 100000
    assert flex.min(a) == 0
    assert flex.max(a) == 9
    a = a.as_double()
    assert approx_equal(flex.mean(a), 4.5, eps=1.e-1)
    assert approx_equal(
      flex.mean(a*a) - flex.mean(a)*flex.mean(a), 8.25, eps=1.e-1)
  for i_trial in xrange(10):
    a = flex.random_double(size=100000, factor=10)
    assert a.size() == 100000
    assert flex.min(a) >= 0
    assert flex.max(a) < 10
    assert approx_equal(flex.mean(a), 5, eps=1.e-1)
    assert approx_equal(
      flex.mean(a*a) - flex.mean(a)*flex.mean(a), 8.25, eps=0.2)
  flex.set_random_seed(value=0)
  assert tuple(flex.random_size_t(3)) \
      == (1791095845, 4282876139L, 3093770124L)
  assert approx_equal(flex.random_double(3),
    (0.9325573593386588, 0.12812444792935673, 0.99904051532414473))
  assert list(flex.random_permutation(size=5)) == [2, 1, 4, 0, 3]
  assert list(flex.random_permutation(size=5)) == [1, 0, 4, 3, 2]
  assert list(flex.random_permutation(size=5)) == [0, 2, 3, 1, 4]
  assert list(flex.random_permutation(size=5)) == [3, 1, 0, 4, 2]
  #
  state = flex.random_generator.getstate()
  r1 = flex.random_size_t(13)
  for i_trial in xrange(10):
    flex.random_generator.setstate(state=state)
    r2 = flex.random_size_t(13)
    assert r2.all_eq(r1)

def exercise_flex_vec3_double():
  flex.exercise_triple(flex.vec3_double, as_double=True)
  a = flex.vec3_double(((1,2,5), (-2,3,4), (3,4,3)))
  assert approx_equal(a.min(), (-2.0,2.0,3.0))
  assert approx_equal(a.max(), (3.0,4.0,5.0))
  assert approx_equal(a.sum(), (2.0,9.0,12.0))
  assert approx_equal(a.mean(), (2.0/3,9.0/3,12.0/3))
  weights = flex.double([1,1,1])
  assert approx_equal(a.mean_weighted(weights=weights), (2.0/3,9.0/3,12.0/3))
  weights = flex.double([2,3,5])
  assert approx_equal(a.mean_weighted(weights=weights), (1.1,3.3,3.7))
  a += (10,20,30)
  assert approx_equal(tuple(a), ((11,22,35), (8,23,34), (13,24,33)))
  assert approx_equal(tuple(a+(20,30,10)), ((31,52,45),(28,53,44),(33,54,43)))
  assert approx_equal(tuple(a+a),
    ((2*11,2*22,2*35), (2*8,2*23,2*34), (2*13,2*24,2*33)))
  a -= (10,20,30)
  assert approx_equal(tuple(a), ((1,2,5), (-2,3,4), (3,4,3)))
  b = a.deep_copy()
  b *= 3
  assert approx_equal(tuple(b), ((3,6,15), (-6,9,12), (9,12,9)))
  b += a
  assert approx_equal(tuple(b), ((4,8,20), (-8,12,16), (12,16,12)))
  assert approx_equal(tuple(a-(20,30,10)),
    ((-19,-28,-5),(-22,-27,-6),(-17,-26,-7)))
  assert tuple(a-a) == ((0,0,0),(0,0,0),(0,0,0))
  a += (10,20,30)
  assert approx_equal(tuple(a*2), ((22,44,70), (16,46,68), (26,48,66)))
  assert approx_equal(tuple(-3*a),
    ((-33,-66,-105), (-24,-69,-102), (-39,-72,-99)))
  assert approx_equal(a/flex.double([1,-2,3]),
    [(11,22,35), (-4,-11.5,-17), (4+1/3.,8,11)])
  assert approx_equal(tuple(a*(-1,1,0,1,0,-1,1,-1,1)),
    ((46,-24,13),(49,-26,11),(44,-20,9)))
  assert approx_equal(tuple((-1,1,0,1,0,-1,1,-1,1)*a),
    ((11,-24,24),(15,-26,19),(11,-20,22)))
  x = flex.double([1,2,3])
  y = flex.double([4,5,6])
  z = flex.double([7,8,9])
  a = flex.vec3_double(x,y,z)
  assert approx_equal(tuple(a), ((1,4,7), (2,5,8), (3,6,9)))
  assert approx_equal(tuple(a.dot(a)), (66,93,126))
  assert approx_equal(tuple(a.dot()), (66,93,126))
  b = flex.vec3_double(z,x,y)
  assert approx_equal(
    a.max_distance(b)**2,
    max(flex.vec3_double(a.as_double()-b.as_double()).dot()))
  assert approx_equal(a.sum_sq(), 285)
  assert approx_equal(a.sum_sq(), flex.sum_sq(a.as_double()))
  assert approx_equal(a.norm(), math.sqrt(285))
  assert approx_equal(a.rms_difference(b), math.sqrt(flex.mean((a-b).dot())))
  assert approx_equal(a.rms_difference(b), b.rms_difference(a))
  assert approx_equal(a.rms_difference(a), 0)
  assert approx_equal(b.rms_difference(b), 0)
  assert approx_equal(a.rms_length(), math.sqrt(flex.mean(a.dot())))
  assert approx_equal((a-a).rms_length(), 0)
  for i_trial in xrange(10):
    for n in [7,10,13]:
      sites_1 = flex.vec3_double(flex.random_double(n*3)*5)
      sites_2 = flex.vec3_double(flex.random_double(n*3)*7)
      m1 = matrix.rec(sites_1.as_double(), (sites_1.size(), 3))
      m2 = matrix.rec(sites_2.as_double(), (sites_2.size(), 3))
      assert approx_equal(
        sites_1.transpose_multiply(sites_2),
        m1.transpose()*m2)
  #
  a = flex.vec3_double([(1,2,5), (-2,3,4), (3,4,3)])
  i = flex.size_t([0,2])
  v = flex.vec3_double([(6,2,-8), (-4,9,2)])
  assert a.add_selected(indices=i, values=v) is a
  assert approx_equal(a, [(7,4,-3), (-2,3,4), (-1,13,5)])

def exercise_histogram():
  x = flex.double(xrange(20))
  h = flex.histogram(data=x)
  assert h.slots().size() == 1000
  h = flex.histogram(x, n_slots=5)
  assert approx_equal(h.data_min(), 0)
  assert approx_equal(h.data_max(), 19)
  assert approx_equal(h.slot_width(), 19/5.)
  assert tuple(h.slots()) == (4,4,4,4,4)
  assert h.n_out_of_slot_range() == 0
  assert approx_equal(h.get_cutoff(max_points=15), 7.60038)
  assert approx_equal(h.get_cutoff(15, relative_tolerance=0.1), 7.98)
  y = flex.double(xrange(-3,23))
  hy = flex.histogram(other=h, data=y)
  assert approx_equal(hy.data_min(), 0)
  assert approx_equal(hy.data_max(), 19)
  assert approx_equal(hy.slot_width(), 19/5.)
  assert tuple(hy.slots()) == (4,4,4,4,4)
  assert hy.n_out_of_slot_range() == 6
  hy = flex.histogram(other=h, data=y, relative_tolerance=0.5)
  assert tuple(hy.slots()) == (5,4,4,4,5)
  assert hy.n_out_of_slot_range() == 4
  hy = flex.histogram(other=h, data=y, relative_tolerance=1)
  assert tuple(hy.slots()) == (7,4,4,4,7)
  assert hy.n_out_of_slot_range() == 0
  s = StringIO()
  hy.show(f=s, prefix="*")
  assert s.getvalue() == """\
*0 - 3.8: 7
*3.8 - 7.6: 4
*7.6 - 11.4: 4
*11.4 - 15.2: 4
*15.2 - 19: 7
"""

def simple_linear_regression(x_obs, y_obs, w_obs):
  assert len(x_obs) == len(y_obs)
  assert len(x_obs) == len(w_obs)
  w = 0
  x = 0
  y = 0
  x2 = 0
  y2 = 0
  xy = 0
  for i in xrange(len(x_obs)):
    w += w_obs[i]
    x += x_obs[i]*w_obs[i]
    y += y_obs[i]*w_obs[i]
    x2 += x_obs[i]**2*w_obs[i]
    y2 += y_obs[i]**2*w_obs[i]
    xy += x_obs[i]*y_obs[i]*w_obs[i]
  determinant = w*x2 - x**2
  a = x2*y - x*xy
  a /= determinant
  b = w*xy - x*y
  b /= determinant
  return a, b

def exercise_linear_regression():
  x = flex.double((1,2,3))
  r = flex.linear_regression(x=x, y=x, epsilon=1.e-6)
  assert r.is_well_defined()
  assert approx_equal(r.y_intercept(), 0)
  assert approx_equal(r.slope(), 1)
  y = flex.double((-1./2+1,-2./2+1,-3./2+1))
  r = flex.linear_regression(x, y)
  assert r.is_well_defined()
  assert approx_equal(r.y_intercept(), 1)
  assert approx_equal(r.slope(), -1./2)
  y = flex.double((0,0,0))
  r = flex.linear_regression(x, y)
  assert r.is_well_defined()
  assert approx_equal(r.y_intercept(), 0)
  assert approx_equal(r.slope(), 0)
  r = flex.linear_regression(y, y)
  assert not r.is_well_defined()
  s = StringIO()
  r.show_summary(f=s)
  assert s.getvalue() == """\
is_well_defined: %s
y_intercept: 0.0
slope: 0.0
""" % str(False)
  y = flex.double((-1./2+1,-2./2+1,-3./2+1))
  for weight in [1,math.pi]:
    weights = flex.double(3, weight)
    r = flex.linear_regression(x=x, y=y, epsilon=1.e-6)
    rw = flex.linear_regression(x=x, y=y, weights=weights, epsilon=1.e-6)
    assert rw.is_well_defined()
    assert approx_equal(rw.y_intercept(), r.y_intercept())
    assert approx_equal(rw.slope(), r.slope())
    a, b = simple_linear_regression(x, y, weights)
    assert approx_equal(a, r.y_intercept())
    assert approx_equal(b, r.slope())
  for i_trial in xrange(100):
    x = flex.random_double(size=10)
    y = x + flex.random_double(size=10) * 0.2 - 0.1 + 3
    weights = flex.random_double(size=10)
    rw = flex.linear_regression(x=x, y=y, weights=weights)
    a, b = simple_linear_regression(x, y, weights)
    assert approx_equal(a, rw.y_intercept())
    assert approx_equal(b, rw.slope())
    assert approx_equal(rw.y_intercept(), 3, 1)

def exercise_linear_correlation():
  x = flex.double((1,2,3))
  c = flex.linear_correlation(x=x, y=x, epsilon=1.e-6)
  assert c.is_well_defined()
  assert c.n() == 3
  assert approx_equal(c.mean_x(), 2)
  assert approx_equal(c.mean_y(), 2)
  assert approx_equal(c.numerator(), 2)
  assert approx_equal(c.sum_denominator_x(), 2)
  assert approx_equal(c.sum_denominator_y(), 2)
  assert approx_equal(c.denominator(), 2)
  assert approx_equal(c.coefficient(), 1)
  y = flex.double((-1./2+1,-2./2+1,-3./2+1))
  c = flex.linear_correlation(x, y)
  assert c.is_well_defined()
  assert approx_equal(c.coefficient(), -1)
  y = flex.double((0,0,0))
  c = flex.linear_correlation(x, y)
  assert c.is_well_defined()
  assert approx_equal(c.coefficient(), 1)
  c = flex.linear_correlation(y, y)
  assert c.is_well_defined()
  c = flex.linear_correlation(flex.double(), flex.double())
  assert not c.is_well_defined()
  s = StringIO()
  c.show_summary(f=s)
  assert s.getvalue() == """\
is_well_defined: %s
mean_x: 0.0
mean_y: 0.0
coefficient: 0.0
""" % str(False)

def exercise_mean_and_variance():
  x = flex.double((1,2,3))
  for w in (None, flex.double((1,1,1))):
    if (w is None):
      mv = flex.mean_and_variance(x)
      assert not mv.have_weights()
    else:
      mv = flex.mean_and_variance(x, w)
      assert mv.have_weights()
    assert approx_equal(mv.mean(), 2)
    assert approx_equal(mv.sum_weights(), 3)
    assert approx_equal(mv.sum_weights_sq(), 3)
    assert approx_equal(mv.sum_weights_values(), 6)
    assert approx_equal(mv.sum_weights_delta_sq(), 2)
    assert approx_equal(mv.gsl_stats_wvariance(), 1)
    assert approx_equal(mv.gsl_stats_wsd(), 1)
    assert approx_equal(
      mv.standard_error_of_mean_calculated_from_sample_weights(),
      1/3**0.5)
    if (not mv.have_weights()):
      assert approx_equal(mv.unweighted_sample_variance(), 1)
      assert approx_equal(mv.unweighted_sample_standard_deviation(), 1)
      assert approx_equal(mv.unweighted_standard_error_of_mean(), 1/3**0.5)
  w = flex.double((1,3,2))
  mv = flex.mean_and_variance(x, w)
  assert approx_equal(mv.mean(), 13/6.)
  assert approx_equal(mv.sum_weights(), 6)
  assert approx_equal(mv.sum_weights_sq(), 14)
  assert approx_equal(mv.sum_weights_values(), 13)
  assert approx_equal(mv.sum_weights_delta_sq(), 17/6.)
  assert approx_equal(mv.gsl_stats_wvariance(),
    (6./(36-14))*((1*(1-13/6.)**2)+(3*(2-13/6.)**2)+(2*(3-13/6.)**2)))
  assert approx_equal(mv.gsl_stats_wsd(), math.sqrt(17/22.))
  assert approx_equal(
    mv.standard_error_of_mean_calculated_from_sample_weights(),
    1/6**0.5)
  x = flex.double((1,1,1))
  mv = flex.mean_and_variance(x)
  assert approx_equal(mv.mean(), 1)
  assert approx_equal(mv.gsl_stats_wvariance(), 0)
  assert approx_equal(
    mv.standard_error_of_mean_calculated_from_sample_weights(),
    1/3**0.5)
  assert approx_equal(mv.unweighted_sample_variance(), 0)
  assert approx_equal(mv.unweighted_sample_standard_deviation(), 0)
  assert approx_equal(mv.unweighted_standard_error_of_mean(), 0)

def exercise_linear_interpolation():
  for flex_type in (flex.float, flex.double):
    tab_x = flex_type([1,2,3,4])
    tab_y = flex_type([10,20,30,40])
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 3.5, 1.e-6),35)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 2.1), 21)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 1), 10)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 1-1.e-6), 10,
      eps=1.e-4)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 1+1.e-6), 10,
      eps=1.e-4)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 4), 40)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 4+1.e-6), 40,
      eps=1.e-4)
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, 4-1.e-6), 40,
      eps=1.e-4)
    for i in xrange(11):
      x = 1+3*i/10.
      assert approx_equal(flex.linear_interpolation(tab_x, tab_y, x), 10*x)
    for x in [0,5]:
      try: flex.linear_interpolation(tab_x, tab_y, 0)
      except KeyboardInterrupt: raise
      except: pass
      else: raise RuntimeError("Exception expected.")
    x = flex_type([1.3,2.4,3.6])
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, x, 1.e-6),
      [13,24,36])
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, x),
      [13,24,36])
    tab_x = flex_type([2,4,8,16])
    tab_y = flex_type([4,6,8,10])
    x = flex_type([2,3,5,6,8,11,12,16])
    assert approx_equal(flex.linear_interpolation(tab_x, tab_y, x),
      [4, 5, 6.5, 7, 8, 8.75, 9, 10])
    tab_x = flex_type([2,4,8,16])
    tab_y = flex_type([4,-6,8,-10])
    for eps in [1.e-6,-1.e-6]:
      x = flex_type([2,3,5,6,8,11,12,16]) + eps
      assert approx_equal(flex.linear_interpolation(tab_x, tab_y, x),
        [4, -1, -6+1/4.*14, -6+2/4.*14, 8, 8-18*3./8, 8-18*4./8, -10],
      eps=1.e-4)

def exercise_loops():
  points = []
  loop = flex.nested_loop(end=(2,3))
  while (not loop.over()):
    points.append(loop())
    loop.incr()
  assert points == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
  points = []
  for index in flex.nested_loop((2,3)):
    points.append(index)
  assert points == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
  points = []
  for index in flex.nested_loop(end=(2,3), open_range=True):
    points.append(index)
  assert points == [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
  points = []
  for index in flex.nested_loop((2,3), False):
    points.append(index)
  assert points == [(0,0),(0,1),(0,2),(0,3),
                    (1,0),(1,1),(1,2),(1,3),
                    (2,0),(2,1),(2,2),(2,3)]
  points = []
  for index in flex.nested_loop(begin=(1,2),end=(3,4)):
    points.append(index)
  assert points == [(1, 2), (1, 3), (2, 2), (2, 3)]
  points = []
  for index in flex.nested_loop(begin=(1,2),end=(3,4),open_range=True):
    points.append(index)
  assert points == [(1, 2), (1, 3), (2, 2), (2, 3)]
  points = []
  for index in flex.nested_loop((1,2),(3,4),False):
    points.append(index)
  assert points == [(1,2),(1,3),(1,4),(2,2),(2,3),(2,4),(3,2),(3,3),(3,4)]
  points = []
  for index in flex.nested_loop((1,2),(1,2),False):
    points.append(index)
  assert points == [(1, 2)]
  points = []
  for index in flex.nested_loop((1,2),(1,2)):
    points.append(index)
  assert points == []
  points = []
  for index in flex.nested_loop([],[]):
    points.append(index)
  assert points == []

def exercise_extract_attributes():
  class group(object):
    def __init__(self, a, b):
      self.a = a
      self.b = b
  groups = []
  for i in xrange(3):
    groups.append(group(i, i+1))
  for array in [groups, tuple(groups)]:
    as = flex.extract_double_attributes(
      array=array, attribute_name="a", none_substitute=None)
    assert approx_equal(as, [0,1,2])
    bs = flex.extract_double_attributes(
      array=array, attribute_name="b", none_substitute=None)
    assert approx_equal(bs, [1,2,3])
  groups[1].a = None
  groups[2].b = None
  for array in [groups, tuple(groups)]:
    as = flex.extract_double_attributes(
      array=array, attribute_name="a", none_substitute=3)
    assert approx_equal(as, [0,3,2])
    bs = flex.extract_double_attributes(
      array=array, attribute_name="b", none_substitute=-4)
    assert approx_equal(bs, [1,2,-4])

def exercise_exceptions():
  f = flex.double(flex.grid((2,3)))
  try: f.assign(1, 0)
  except RuntimeError, e:
    assert str(e) == "Array must be 0-based 1-dimensional."
  else:
    raise AssertionError, "No exception or wrong exception."
  try: f.push_back(0)
  except RuntimeError, e:
    assert str(e) == "Array must be 0-based 1-dimensional."
  else:
    raise AssertionError, "No exception or wrong exception."
  try: f[(2,0)]
  except IndexError, e:
    assert str(e) == "Index out of range."
  else:
    raise AssertionError, "No exception or wrong exception."

def exercise_matrix():
  for ag,bg in [[(0,0),(0,0)],[(0,1),(1,0)],[(1,0),(0,2)]]:
    a = flex.double(flex.grid(ag))
    b = flex.double(flex.grid(bg))
    c = a.matrix_multiply(b)
    assert c.focus() == (ag[0],bg[1])
    assert c.all_eq(0)
    for m in [a,b]:
      mtm = m.matrix_transpose_multiply_as_packed_u()
      n = m.focus()[1]
      assert mtm.size() == n*(n+1)/2
      assert mtm.all_eq(0)
  a = flex.double(range(1,7))
  a.resize(flex.grid(3,2))
  b = flex.double(range(1,7))
  b.resize(flex.grid(2,3))
  c = a.matrix_multiply(b)
  assert c.focus() == (3,3)
  assert list(c) == [9, 12, 15, 19, 26, 33, 29, 40, 51]
  for a_n_rows in xrange(1,4):
    for a_n_columns in xrange(1,4):
      a = flex.random_double(size=a_n_rows*a_n_columns)
      b = flex.random_double(size=a_n_rows*a_n_columns)
      c = a.matrix_multiply(b)
      d = matrix.row(a) * matrix.col(b)
      assert approx_equal(c, d[0])
      assert approx_equal(a.dot(b), matrix.row(a).dot(matrix.col(b)))
      a.resize(flex.grid(a_n_rows, a_n_columns))
      for b_n_columns in xrange(1,4):
        b = flex.random_double(size=a_n_columns*b_n_columns)
        b.resize(flex.grid(a_n_columns, b_n_columns))
        c = a.matrix_multiply(b)
        d = matrix.rec(a, a.focus()) * matrix.rec(b, b.focus())
        assert c.focus() == d.n
        assert approx_equal(c, d)
        ata = a.matrix_transpose_multiply_as_packed_u() \
          .matrix_packed_u_as_symmetric()
        assert approx_equal(ata, a.matrix_transpose().matrix_multiply(a))
        assert c.matrix_transpose().focus() == (d.n[1], d.n[0])
        assert approx_equal(c.matrix_transpose(), d.transpose())
        assert approx_equal(c.matrix_transpose().matrix_transpose(), d)
        c.matrix_transpose_in_place()
        assert c.focus() == (d.n[1], d.n[0])
        assert approx_equal(c, d.transpose())
        c.matrix_transpose_in_place()
        assert c.focus() == d.n
        assert approx_equal(c, d)
        b = flex.random_double(size=a_n_columns)
        c = a.matrix_multiply(b)
        d = matrix.rec(a, a.focus()) * matrix.col(b)
        assert d.n[1] == 1
        assert c.focus() == (d.n[0],)
        assert approx_equal(c, d)
        b = flex.random_double(size=a_n_rows)
        c = b.matrix_multiply(a)
        d = matrix.row(b) * matrix.rec(a, a.focus())
        assert d.n[0] == 1
        assert c.focus() == (d.n[1],)
        assert approx_equal(c, d)
  #
  for n_rows in xrange(10):
    for n_columns in xrange(10):
      m = flex.random_double(size=n_rows*n_columns)*4-2
      m.reshape(flex.grid(n_rows,n_columns))
      p = flex.random_double(size=n_columns*(n_columns+1)/2)
      nopt = m.matrix_multiply(p.matrix_packed_u_as_symmetric())
      assert nopt.focus() == (n_rows, n_columns)
      opt = m.matrix_multiply_packed_u(p)
      assert approx_equal(opt, nopt)
      nopt = nopt.matrix_multiply(m.matrix_transpose()) \
        .matrix_symmetric_as_packed_u()
      assert nopt.size() == n_rows*(n_rows+1)/2
      opt = m.matrix_multiply_packed_u_multiply_lhs_transpose(packed_u=p)
      assert approx_equal(opt, nopt)
      #
      p = flex.complex_double(
        reals=flex.random_double(size=n_columns*(n_columns+1)/2),
        imags=flex.random_double(size=n_columns*(n_columns+1)/2))
      ps = p.matrix_packed_u_as_symmetric()
      nopt = m.matrix_multiply(ps)
      assert nopt.focus() == (n_rows, n_columns)
      opt = m.matrix_multiply_packed_u(p)
      assert approx_equal(opt, nopt)
      nopts = nopt.matrix_multiply(m.matrix_transpose())
      nopt = flex.complex_double(
        reals=flex.real(nopts).matrix_symmetric_as_packed_u(),
        imags=flex.imag(nopts).matrix_symmetric_as_packed_u())
      assert nopt.size() == n_rows*(n_rows+1)/2
      opt = m.matrix_multiply_packed_u_multiply_lhs_transpose(packed_u=p)
      assert approx_equal(opt, nopt)
  #
  a = flex.polar(flex.double(range(1,6+1)), flex.double(range(2,7+1)))
  a.resize(flex.grid(2,3))
  b = a.deep_copy()
  b.resize(flex.grid(3,2))
  c = a.matrix_multiply(b)
  assert approx_equal(c, matrix.rec(a, a.focus())*matrix.rec(b, b.focus()))
  b = flex.real(b)
  c = a.matrix_multiply(b)
  assert approx_equal(c, matrix.rec(a, a.focus())*matrix.rec(b, b.focus()))
  c = b.matrix_multiply(a)
  assert approx_equal(c, matrix.rec(b, b.focus())*matrix.rec(a, a.focus()))
  #
  assert flex.double().matrix_outer_product(rhs=flex.double()).size() == 0
  assert flex.double([1]).matrix_outer_product(rhs=flex.double()).size() == 0
  assert flex.double().matrix_outer_product(rhs=flex.double([1])).size() == 0
  op = flex.double([1,2,3]).matrix_outer_product(rhs=flex.double([4,5]))
  assert op.focus() == (3,2)
  assert approx_equal(op, [4,5,8,10,12,15])
  lu = flex.double([1,0,0,0,1,0,0,0,1])
  lu.resize(flex.grid(3,3))
  pivot_indices = lu.matrix_lu_decomposition_in_place()
  assert approx_equal(lu, [1,0,0,0,1,0,0,0,1])
  assert list(pivot_indices) == [0,1,2,0]
  for rank in range(1,6) + [18]:
    for i in xrange(10):
      m = flex.random_double(size=rank*rank)*3-1.5
      m.resize(flex.grid(rank,rank))
      assert approx_equal(
        m.matrix_diagonal(), [m[j*rank+j] for j in xrange(rank)])
      mc = m.deep_copy()
      mc.matrix_diagonal_set_in_place(value=3)
      assert approx_equal(mc.matrix_diagonal(), [3]*rank)
      mc.matrix_diagonal_add_in_place(value=-5)
      assert approx_equal(mc.matrix_diagonal(), [-2]*rank)
      assert approx_equal(
        m.matrix_diagonal_sum(), flex.sum(m.matrix_diagonal()))
      assert m.matrix_trace() == m.matrix_diagonal_sum()
      assert approx_equal(
        m.matrix_diagonal_product(), flex.product(m.matrix_diagonal()))
      b = flex.random_double(size=rank)*7-3.5
      lu = m.deep_copy()
      pivot_indices = lu.matrix_lu_decomposition_in_place()
      x = lu.matrix_lu_back_substitution(pivot_indices=pivot_indices, b=b)
      assert approx_equal(m.matrix_multiply(x), b)
  lu = flex.double([0,0,0,0,1,0,0,0,1])
  lu.resize(flex.grid(3,3))
  try: lu.matrix_lu_decomposition_in_place()
  except RuntimeError, e:
    assert str(e) == "lu_decomposition_in_place: singular matrix"
  else: raise RuntimeError("Exception expected.")
  lu = flex.double([1,0,0,0,1,0,0,0,1])
  lu.resize(flex.grid(3,3))
  b = flex.double([1,2,3])
  pivot_indices = flex.size_t([0,1,4,0])
  try: lu.matrix_lu_back_substitution(pivot_indices, b)
  except RuntimeError, e:
    assert str(e) == "lu_back_substitution: pivot_indices[i] out of range"
  else: raise RuntimeError("Exception expected.")
  lu = flex.double([1,6,4,32,6,2,-1,63,-4,1,4,6,1,0,-13,5])
  lu.resize(flex.grid(4,4))
  pivot_indices = lu.matrix_lu_decomposition_in_place()
  assert approx_equal(
    lu.matrix_determinant_via_lu(pivot_indices=pivot_indices), -16522)
  m = [1,6,4,32,6,2,-1,63,-4,1,4,6,1,0,-13,5]
  assert approx_equal(matrix.sqr(m).determinant(), -16522)
  m = [1,6,4,32,6,2,-1,-63,-4,1,4,6,1,0,-13,5]
  assert approx_equal(matrix.sqr(m).determinant(), 21908)
  m = [0,0,0,0,6,2,-1,63,-4,1,4,6,1,0,-13,5]
  assert matrix.sqr(m).determinant() == 0
  #
  for size in xrange(10):
    a = flex.double(size)
    b = flex.double(size)
    assert a.cos_angle(b=b, value_if_undefined=-10) == -10
    assert a.cos_angle(b=b) is None
    assert a.angle(b=b) is None
  for size in xrange(1,10):
    a = flex.double(xrange(10))
    b = flex.double(xrange(10))
    assert approx_equal(a.cos_angle(b=b, value_if_undefined=-10), 1,eps=1.e-10)
    assert approx_equal(a.cos_angle(b=b), 1, eps=1.e-10)
    assert approx_equal(a.angle(b=b), 0, eps=1.e-10)
  a = flex.double([1,0])
  b = flex.double([1,1])
  assert approx_equal(a.angle(b, deg=True), 45)
  b = flex.double([1,2])
  assert approx_equal(a.angle(b, deg=True), 63.4349488229)
  mersenne_twister = flex.mersenne_twister(0)
  for size in xrange(1,100):
    a = mersenne_twister.random_double(size=size)
    b = mersenne_twister.random_double(size=size)
    assert approx_equal(
      a.cos_angle(b=b), matrix.col(a).cos_angle(matrix.col(b)))
    assert approx_equal(a.angle(b=b, deg=False), a.angle(b))
    assert approx_equal(a.angle(b=b, deg=True), a.angle(b)*180/math.pi)
  #
  assert approx_equal(
    flex.double([[1,2,3],[4,5,6],[7,8,9]]).matrix_upper_triangle_as_packed_u(),
    [1,2,3,5,6,9])
  assert approx_equal(
    flex.double([1,2,3,5,6,9]).matrix_packed_u_as_upper_triangle(),
    [1,2,3, 0,5,6, 0,0,9])
  assert approx_equal(
    flex.double([[1,2,3],[4,5,6],[7,8,9]]).matrix_lower_triangle_as_packed_l(),
    [1,4,5,7,8,9])
  assert approx_equal(
    flex.double([1,4,5,7,8,9]).matrix_packed_l_as_lower_triangle(),
    [1,0,0, 4,5,0, 7,8,9])
  m = flex.double(xrange(1,17))
  m.resize(flex.grid(4, 4))
  assert approx_equal(
    m.matrix_upper_triangle_as_packed_u(),
    [1,2,3,4,6,7,8,11,12,16])
  assert approx_equal(
    flex.double([1,2,3,4,6,7,8,11,12,16])
      .matrix_packed_u_as_upper_triangle(),
    [1,2,3,4,0,6,7,8,0,0,11,12,0,0,0,16])
  assert approx_equal(
    m.matrix_lower_triangle_as_packed_l(),
    [1,5,6,9,10,11,13,14,15,16])
  assert approx_equal(
    flex.double([1,5,6,9,10,11,13,14,15,16])
      .matrix_packed_l_as_lower_triangle(),
    [1,0,0,0,5,6,0,0,9,10,11,0,13,14,15,16])
  def exercise_packed(n, s, u, l):
    assert s.focus() == (n,n)
    if (u is not None):
      p = s.matrix_upper_triangle_as_packed_u()
      assert approx_equal(p, u)
      t = p.matrix_packed_u_as_upper_triangle()
      assert approx_equal(t.matrix_upper_triangle_as_packed_u(), u)
      p = s.matrix_symmetric_as_packed_u()
      assert approx_equal(p, u)
      us = p.matrix_packed_u_as_symmetric()
      assert us.focus() == (n,n)
      assert approx_equal(us, s)
    if (l is not None):
      p = s.matrix_lower_triangle_as_packed_l()
      assert approx_equal(p, l)
      t = p.matrix_packed_l_as_lower_triangle()
      assert approx_equal(t.matrix_lower_triangle_as_packed_l(), l)
      p = s.matrix_symmetric_as_packed_l()
      assert approx_equal(p, l)
      ls = p.matrix_packed_l_as_symmetric()
      assert ls.focus() == (n,n)
      assert approx_equal(ls, s)
  exercise_packed(
    n=0,
    s=flex.double(flex.grid(0,0)),
    u=flex.double(),
    l=flex.double())
  exercise_packed(
    n=1,
    s=flex.double([[1]]),
    u=flex.double([1]),
    l=flex.double([1]))
  exercise_packed(
    n=2,
    s=flex.double([[1,2],[2,3]]),
    u=flex.double([1,2,3]),
    l=flex.double([1,2,3]))
  exercise_packed(
    n=3,
    s=flex.double([[1,2,3],[2,4,5],[3,5,6]]),
    u=flex.double([1,2,3,4,5,6]),
    l=flex.double([1,2,4,3,5,6]))
  exercise_packed(
    n=4,
    s=flex.double([[1,-2,3,4], [-2,-5,6,-7], [3,6,8,9], [4,-7,9,-10]]),
    u=flex.double([1,-2,3,4,-5,6,-7,8,9,-10]),
    l=flex.double([1,-2,-5,3,6,8,4,-7,9,-10]))
  for n in xrange(20):
    p = flex.random_double(size=n*(n+1)/2)
    exercise_packed(
      n=n,
      s=p.matrix_packed_u_as_symmetric(),
      u=p,
      l=None)
    exercise_packed(
      n=n,
      s=p.matrix_packed_l_as_symmetric(),
      u=None,
      l=p)
  e = flex.double([[1,2,3],[2,4,5],[3,9,6]])
  try: e.matrix_symmetric_as_packed_u()
  except RuntimeError, err:
    assert str(err) == "symmetric_as_packed_u(): matrix is not symmetric."
  else: raise RuntimeError("Exception expected.")
  p = e.matrix_symmetric_as_packed_u(relative_epsilon=1)
  assert approx_equal(p, [1,2,3,4,7,6])
  try: e.matrix_symmetric_as_packed_l()
  except RuntimeError, err:
    assert str(err) == "symmetric_as_packed_l(): matrix is not symmetric."
  else: raise RuntimeError("Exception expected.")
  p = e.matrix_symmetric_as_packed_l(relative_epsilon=1)
  assert approx_equal(p, [1,2,4,3,7,6])
  p = flex.double(4)
  for method in [p.matrix_packed_u_as_symmetric,
                 p.matrix_packed_l_as_symmetric]:
    try: method()
    except RuntimeError, e:
      assert str(e).endswith(
        "SCITBX_ASSERT(n*(n+1)/2 == packed_size) failure.")
    else: raise RuntimeError("Exception expected.")
  #
  from scitbx.examples import immoptibox_ports
  immoptibox_ports.py_cholesky_decomposition \
    = immoptibox_ports.cholesky_decomposition
  immoptibox_ports.cholesky_decomposition = exercise_cholesky_decomposition
  immoptibox_ports.tst_flex_counts = 0
  immoptibox_ports.exercise_cholesky()
  immoptibox_ports.cholesky_decomposition \
    = immoptibox_ports.py_cholesky_decomposition
  assert immoptibox_ports.tst_flex_counts == 299
  del immoptibox_ports.tst_flex_counts

def exercise_cholesky_decomposition(a):
  from scitbx.examples import immoptibox_ports
  c = immoptibox_ports.py_cholesky_decomposition(a)
  au = a.matrix_symmetric_as_packed_u()
  cl = au.matrix_cholesky_decomposition()
  if (c is None):
    assert cl.size() == 0
  else:
    assert approx_equal(cl, c.matrix_lower_triangle_as_packed_l())
    cu = cl.matrix_packed_l_as_symmetric().matrix_symmetric_as_packed_u()
    for i_trial in xrange(10):
      b = flex.random_double(size=a.focus()[0], factor=2)-1
      x = cu.matrix_cholesky_solve_packed_u(b=b)
      assert approx_equal(a.matrix_multiply(x), b)
      if (i_trial == 0):
        pivots = flex.size_t(xrange(b.size()))
      else:
        pivots = flex.random_permutation(size=b.size())
      xp = cu.matrix_cholesky_solve_packed_u(
        b=b.select(pivots, reverse=True), pivots=pivots).select(pivots)
      assert approx_equal(xp, x)
  immoptibox_ports.tst_flex_counts += 1
  return c

def exercise_matrix_cholesky_gill_murray_wright():
  import scitbx.math
  def p_as_mx(p):
    n = len(p)
    m = [0]*n**2
    for i in xrange(n):
      m[p[i]*n+i] = 1
    return matrix.sqr(m)
  def core(a):
    c = flex.double(a)
    c.resize(flex.grid(a.n))
    u = c.matrix_upper_triangle_as_packed_u()
    gwm = u.matrix_cholesky_gill_murray_wright_decomposition_in_place(
      epsilon=1.e-8)
    assert gwm.epsilon == 1.e-8
    u = c.matrix_upper_triangle_as_packed_u()
    gwm = u.matrix_cholesky_gill_murray_wright_decomposition_in_place()
    assert gwm.epsilon == scitbx.math.floating_point_epsilon_double_get()
    assert gwm.packed_u.id() == u.id()
    p, e = gwm.pivots, gwm.e
    r = matrix.sqr(u.matrix_packed_u_as_upper_triangle())
    rtr = r.transpose() * r
    pm = p_as_mx(p)
    ptap = pm.transpose() * a * pm
    ptaep = ptap + matrix.diag(e)
    assert approx_equal(ptaep, rtr)
    b = flex.random_double(size=a.n[0], factor=2)-1
    x = gwm.solve(b=b)
    ae = pm * ptaep * pm.transpose()
    assert approx_equal(ae*matrix.col(x), b)
    return p, e, r
  # empty matrix
  a = matrix.sqr([])
  p, e, r = core(a)
  assert p.size() == 0
  assert e.size() == 0
  assert len(r) == 0
  n_max = 15
  n_trials_per_n = 10
  # identity matrices
  for n in xrange(1,n_max+1):
    a = matrix.diag([1]*n)
    p, e, r = core(a)
    assert list(p) == range(n)
    assert approx_equal(e, [0]*n)
    assert approx_equal(r, a)
  # null matrices
  for n in xrange(1,n_max+1):
    a = matrix.sqr([0]*n*n)
    p, e, r = core(a)
    assert list(p) == range(n)
    assert list(e) == [scitbx.math.floating_point_epsilon_double_get()]*n
    for i in xrange(n):
      for j in xrange(n):
        if (i != j): r(i,j) == 0
        else: r(i,j) == r(0,0)
  # random semi-positive diagonal matrices
  for n in xrange(1,n_max+1):
    for i_trial in xrange(n_trials_per_n):
      a = matrix.diag(flex.random_double(size=n))
      p, e, r = core(a)
      assert approx_equal(e, [0]*n)
      for i in xrange(n):
        for j in xrange(n):
          if (i != j): approx_equal(r(i,j), 0)
  # random diagonal matrices
  for n in xrange(1,n_max+1):
    for i_trial in xrange(n_trials_per_n):
      a = matrix.diag(flex.random_double(size=n, factor=2)-1)
      p, e, r = core(a)
      for i in xrange(n):
        for j in xrange(n):
          if (i != j): approx_equal(r(i,j), 0)
  # random semi-positive definite matrices
  for n in xrange(1,n_max+1):
    for i_trial in xrange(n_trials_per_n):
      m = matrix.sqr(flex.random_double(size=n*n, factor=2)-1)
      a = m.transpose_multiply()
      p, e, r = core(a)
      assert approx_equal(e, [0]*n)
  # random matrices
  for n in xrange(1,n_max+1):
    size = n*(n+1)/2
    for i_trial in xrange(n_trials_per_n):
      a = (flex.random_double(size=size, factor=2)-1) \
            .matrix_packed_u_as_symmetric()
      core(matrix.sqr(a))
      a.matrix_diagonal_set_in_place(0)
      core(matrix.sqr(a))
  # J. Nocedal and S. Wright:
  # Numerical Optimization.
  # Springer, New York, 1999, pp. 145-150.
  for i in xrange(3):
    for j in xrange(3):
      a = flex.double([[4,2,1],[2,6,3],[1,3,-0.004]])
      a.matrix_swap_rows_in_place(i=i, j=j)
      a.matrix_swap_columns_in_place(i=i, j=j)
      p, e, r = core(matrix.sqr(a))
      if (i == 0 and j == 0):
        assert list(p) == [1,0,2]
      assert approx_equal(e, [0.0, 0.0, 3.008])
      assert approx_equal(r,
        [2.4494897427831779, 0.81649658092772592, 1.2247448713915889,
         0.0, 1.8257418583505538, 0.0,
         0.0, 0.0, 1.2263767773404712])

def exercise_matrix_move():
  a = flex.double(flex.grid(0,0))
  b = a.matrix_copy_block(i_row=0, i_column=0, n_rows=0, n_columns=0)
  assert b.focus() == (0,0)
  a.matrix_paste_block_in_place(block=b, i_row=0, i_column=0)
  for a in [flex.double([[1]]), flex.double([[1,2],[3,4]])]:
    b = a.matrix_copy_block(i_row=0, i_column=0, n_rows=0, n_columns=0)
    assert b.focus() == (0,0)
    a.matrix_paste_block_in_place(block=b, i_row=0, i_column=0)
    b = a.matrix_copy_block(i_row=0, i_column=0, n_rows=1, n_columns=0)
    assert b.focus() == (1,0)
    a.matrix_paste_block_in_place(block=b, i_row=0, i_column=0)
    b = a.matrix_copy_block(i_row=0, i_column=0, n_rows=0, n_columns=1)
    assert b.focus() == (0,1)
    a.matrix_paste_block_in_place(block=b, i_row=0, i_column=0)
    b = a.matrix_copy_block(i_row=0, i_column=0, n_rows=1, n_columns=1)
    assert b.focus() == (1,1)
    assert list(b) == [1]
    b[0] = 5
    a.matrix_paste_block_in_place(block=b, i_row=0, i_column=0)
    assert a[0] == 5
  b = a.matrix_copy_block(0,0,2,2)
  assert b.focus() == (2,2)
  assert b.all_eq(a)
  a = flex.double(xrange(1,20+1))
  a.resize(flex.grid(4,5))
  b = a.matrix_copy_block(0,0,4,5)
  assert b.focus() == (4,5)
  assert b.all_eq(a)
  for i in xrange(4):
    for j in xrange(5):
      b = a.matrix_copy_block(i,j,1,1)
      assert b.focus() == (1,1)
      assert b[0] == a[(i,j)]
      c = a.deep_copy()
      c.matrix_paste_block_in_place(flex.double([[93]]), i, j)
      assert c[(i,j)] == 93
  for i in xrange(3):
    for j in xrange(5):
      b = a.matrix_copy_block(i,j,2,1)
      assert b.focus() == (2,1)
      assert b[0] == a[(i,j)]
      assert b[1] == a[(i+1,j)]
      c = a.deep_copy()
      c.matrix_paste_block_in_place(flex.double([[91],[-10]]), i, j)
      assert c[(i,j)] == 91
      assert c[(i+1,j)] == -10
  for i in xrange(3):
    for j in xrange(4):
      b = a.matrix_copy_block(i,j,2,2)
      assert b.focus() == (2,2)
      assert b[0] == a[(i,j)]
      assert b[1] == a[(i,j+1)]
      assert b[2] == a[(i+1,j)]
      assert b[3] == a[(i+1,j+1)]
      c = a.deep_copy()
      c.matrix_paste_block_in_place(flex.double([[97,-4],[-13,84]]), i, j)
      assert c[(i,j)] == 97
      assert c[(i,j+1)] == -4
      assert c[(i+1,j)] == -13
      assert c[(i+1,j+1)] == 84
  for i in xrange(3):
    for j in xrange(3):
      b = a.matrix_copy_block(i,j,2,3)
      assert b.focus() == (2,3)
      assert b[0] == a[(i,j)]
      assert b[1] == a[(i,j+1)]
      assert b[2] == a[(i,j+2)]
      assert b[3] == a[(i+1,j)]
      assert b[4] == a[(i+1,j+1)]
      assert b[5] == a[(i+1,j+2)]
      c = a.deep_copy()
      c.matrix_paste_block_in_place(flex.double([[79,-3,75],[-31,48,-7]]),i,j)
      assert c[(i,j)] == 79
      assert c[(i,j+1)] == -3
      assert c[(i,j+2)] == 75
      assert c[(i+1,j)] == -31
      assert c[(i+1,j+1)] == 48
      assert c[(i+1,j+2)] == -7
  a = flex.complex_double([1,2j,3j,4])
  a.resize(flex.grid(2,2))
  b = a.matrix_copy_block(i_row=1, i_column=0, n_rows=1, n_columns=2)
  assert approx_equal(b, [3j,(4+0j)])
  b = flex.complex_double([10,20j])
  b.resize(flex.grid(2,1))
  assert b.matrix_transpose().focus() == (1,2)
  a.matrix_paste_block_in_place(block=b, i_row=0, i_column=1)
  assert approx_equal(a, [(1+0j), (10+0j), 3j, 20j])
  #
  a = flex.double(flex.grid(0,0))
  a.matrix_copy_upper_to_lower_triangle_in_place()
  a.matrix_copy_lower_to_upper_triangle_in_place()
  a = flex.double([[1]])
  a.matrix_swap_rows_in_place(i=0, j=0)
  a.matrix_swap_columns_in_place(i=0, j=0)
  a.matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place(i=0, j=0)
  assert list(a) == [1]
  u = flex.double([1])
  u.matrix_packed_u_swap_rows_and_columns_in_place(i=0, j=0)
  assert list(u) == [1]
  a = flex.double([[1,2],[3,4]])
  a.matrix_swap_rows_in_place(i=0, j=0)
  assert list(a) == [1,2,3,4]
  a.matrix_swap_rows_in_place(i=0, j=1)
  assert list(a) == [3,4,1,2]
  a.matrix_swap_rows_in_place(i=1, j=0)
  assert list(a) == [1,2,3,4]
  a.matrix_swap_columns_in_place(i=0, j=0)
  assert list(a) == [1,2,3,4]
  a.matrix_swap_columns_in_place(i=0, j=1)
  assert list(a) == [2,1,4,3]
  a.matrix_swap_columns_in_place(i=1, j=0)
  assert list(a) == [1,2,3,4]
  a.matrix_swap_rows_in_place(i=1, j=0)
  a.matrix_swap_columns_in_place(i=1, j=0)
  assert list(a) == [4,3,2,1]
  a.matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place(i=0, j=1)
  assert list(a) == [1,3,2,4]
  u = a.matrix_upper_triangle_as_packed_u()
  u.matrix_packed_u_swap_rows_and_columns_in_place(i=0, j=1)
  assert list(u) == [4,3,1]
  a.matrix_copy_upper_to_lower_triangle_in_place()
  assert list(a) == [1,3,3,4]
  a = flex.double([[1,2],[3,4]])
  a.matrix_copy_lower_to_upper_triangle_in_place()
  assert list(a) == [1,3,3,4]
  a = flex.double([[1,2,3],[4,5,6]])
  a.matrix_swap_rows_in_place(i=1, j=0)
  assert list(a) == [4,5,6,1,2,3]
  a.matrix_swap_columns_in_place(i=1, j=2)
  assert list(a) == [4,6,5,1,3,2]
  a.matrix_swap_columns_in_place(i=2, j=0)
  assert list(a) == [5,6,4,2,3,1]
  a.resize(flex.grid(3,2))
  a.matrix_swap_rows_in_place(i=1, j=2)
  assert list(a) == [5,6,3,1,4,2]
  a.matrix_swap_rows_in_place(i=2, j=0)
  assert list(a) == [4,2,3,1,5,6]
  for n in xrange(1,11): # must be at least 8 to reveal all bugs
    a0 = flex.double(xrange(1,n+1))
    a0.resize(flex.grid(n,n))
    for i in xrange(n):
      for j in xrange(n):
        for triangle_flag in ["u","l"]:
          a = a0.deep_copy()
          if (triangle_flag == "u"):
            a.matrix_copy_upper_to_lower_triangle_in_place()
          else:
            a.matrix_copy_lower_to_upper_triangle_in_place()
          a.matrix_swap_rows_in_place(i=i, j=j)
          a.matrix_swap_columns_in_place(i=i, j=j)
          b = a0.deep_copy()
          if (triangle_flag == "l"):
            b.matrix_transpose_in_place()
          b0 = b.deep_copy()
          b.matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place(
            i=i, j=j)
          for ii in xrange(1,n):
            for jj in xrange(i):
              assert b[(ii,jj)] == b0[(ii,jj)]
          b.matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place(
            i=i, j=j)
          assert list(b) == list(b0)
          b.matrix_symmetric_upper_triangle_swap_rows_and_columns_in_place(
            i=i, j=j)
          b.matrix_copy_upper_to_lower_triangle_in_place()
          assert list(b) == list(a)
          u = b0.matrix_upper_triangle_as_packed_u()
          u.matrix_packed_u_swap_rows_and_columns_in_place(i=i, j=j)
          assert list(u) == list(b.matrix_upper_triangle_as_packed_u())
          u.matrix_packed_u_swap_rows_and_columns_in_place(i=i, j=j)
          assert list(u) == list(b0.matrix_upper_triangle_as_packed_u())

def exercise_matrix_inversion_in_place():
  m = flex.double()
  m.resize(flex.grid(0,0))
  m.matrix_inversion_in_place(m)
  b = flex.double()
  b.resize(flex.grid(0,0))
  m.matrix_inversion_in_place(b=b)
  m = flex.double([2])
  m.resize(flex.grid(1,1))
  m.matrix_inversion_in_place()
  assert approx_equal(m, [1/2.])
  m = flex.double([2,0,0,-3])
  m.resize(flex.grid(2,2))
  m.matrix_inversion_in_place()
  assert approx_equal(m, [1/2.,0,0,-1/3.])
  m = flex.double([1,2,-3,4])
  m.resize(flex.grid(2,2))
  m.matrix_inversion_in_place()
  assert approx_equal(m, [2/5.,-1/5.,3/10.,1/10.])
  m = flex.double([2,0,0,0,-3,0,0,0,4])
  m.resize(flex.grid(3,3))
  m.matrix_inversion_in_place()
  assert approx_equal(m, [1/2.,0,0,0,-1/3.,0,0,0,1/4.])
  m = flex.double([1,2,-3,-2,4,-1,8,0,4])
  m.resize(flex.grid(3,3))
  m.matrix_inversion_in_place()
  assert approx_equal(m, [1/7.,-1/14.,5/56.,0,1/4.,1/16.,-2/7.,1/7.,1/14.])
  from scitbx import matrix
  for n in xrange(1,12):
    u = flex.double(n*n, 0)
    for i in xrange(0,n*n,n+1): u[i] = 1
    for diag in [1,2]:
      m = flex.double(n*n, 0)
      for i in xrange(0,n*n,n+1): m[i] = diag
      m.resize(flex.grid(n,n))
      m_orig = matrix.rec(m, (n,n))
      m.matrix_inversion_in_place()
      m_inv = matrix.rec(m, (n,n))
      assert approx_equal(m_orig*m_inv, u)
      assert approx_equal(m_inv*m_orig, u)
      for n_b in xrange(0,4):
        m = flex.double(m_orig)
        m.resize(flex.grid(n,n))
        b = flex.double(xrange(1,n*n_b+1))
        b.resize(flex.grid(n_b,n))
        b_orig = matrix.rec(b, (n,n_b))
        m.matrix_inversion_in_place(b)
        m_inv = matrix.rec(m, (n,n))
        x = matrix.rec(b, (n_b,n))
        assert approx_equal(m_orig*m_inv, u)
        assert approx_equal(m_inv*m_orig, u)
        for i_b in xrange(n_b):
          b_i = matrix.col(b_orig.elems[i_b*n:(i_b+1)*n])
          x_i = matrix.col(x.elems[i_b*n:(i_b+1)*n])
          assert approx_equal(m_orig*x_i, b_i)
  for n in xrange(1,12):
    u = flex.double(n*n, 0)
    for i in xrange(0,n*n,n+1): u[i] = 1
    for i_trial in xrange(3):
      m = 2*flex.random_double(n*n)-1
      m.resize(flex.grid(n,n))
      m_orig = matrix.rec(m, (n,n))
      try:
        m.matrix_inversion_in_place()
      except RuntimeError, e:
        assert str(e) == "inversion_in_place: singular matrix"
      else:
        m_inv = matrix.rec(m, (n,n))
        assert approx_equal(m_orig*m_inv, u)
        assert approx_equal(m_inv*m_orig, u)
        for n_b in xrange(0,4):
          m = flex.double(m_orig)
          m.resize(flex.grid(n,n))
          b = flex.random_double(n*n_b)
          b.resize(flex.grid(n_b,n))
          b_orig = matrix.rec(b, (n,n_b))
          m.matrix_inversion_in_place(b)
          m_inv = matrix.rec(m, (n,n))
          x = matrix.rec(b, (n_b,n))
          assert approx_equal(m_orig*m_inv, u)
          assert approx_equal(m_inv*m_orig, u)
          for i_b in xrange(n_b):
            b_i = matrix.col(b_orig.elems[i_b*n:(i_b+1)*n])
            x_i = matrix.col(x.elems[i_b*n:(i_b+1)*n])
            assert approx_equal(m_orig*x_i, b_i)

def exercise_pickle_single_buffered():
  a = flex.bool((1,0,1))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,0,1)
  a = flex.double(())
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 0
  a = flex.double((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = flex.int((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = flex.float((1,2,3))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1,2,3)
  a = flex.complex_double((1+2j, 2+3j, 4+5j))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (1+2j, 2+3j, 4+5j)
  a = flex.double(flex.grid((-1,2,-3), (7,5,3)).set_focus((1,3,2)), 13)
  a[(1,2,-1)] = -8
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 8 * 3 * 6
  assert tuple(a) == tuple(b)
  assert a.origin() == b.origin()
  assert a.all() == b.all()
  assert a.focus() == b.focus()
  assert a.accessor() == b.accessor()

def exercise_pickle_double_buffered():
  a = flex.std_string()
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 0
  assert tuple(b) == ()
  a = flex.std_string(("abc", "bcd", "cde"))
  p = pickle.dumps(a)
  b = pickle.loads(p)
  assert b.size() == 3
  assert tuple(b) == (("abc", "bcd", "cde"))

def pickle_large_arrays(max_exp, verbose):
  try:
    for array_type in (
        flex.bool,
        flex.size_t,
        flex.int,
        flex.long,
        flex.float,
        flex.double,
        flex.complex_double,
        flex.std_string):
      for e in xrange(max_exp+1):
        n = 2**e
        if (array_type == flex.bool):
          val = 1
        elif (array_type == flex.size_t):
          val = 2147483647
        elif (array_type == flex.int):
          val = -2147483647
        elif (array_type == flex.long):
          val = -9223372036854775808
          if (type(val) == type(1L)):
            val = -2147483647
        elif (array_type == flex.complex_double):
          val = complex(-1.234567809123456e+20, -1.234567809123456e+20)
        elif (array_type in (flex.float, flex.double)):
          val = -1.234567809123456e+20
        elif (array_type == flex.std_string):
          val = "x" * 10
        else:
          raise AssertionError, "Unexpected array type."
        a = array_type(n, val)
        for pickler in (0, pickle, cPickle):
          if (pickler == 0):
            pickler_name = "g/setstate"
            t0 = time.time()
            s = a.__getstate__()
            td = time.time() - t0
            t0 = time.time()
            b = array_type()
            b.__setstate__(s)
            tl = time.time() - t0
          else:
            pickler_name = pickler.__name__
            f = open("pickle.tmp", "wb")
            t0 = time.time()
            pickler.dump(a, f, 1)
            td = time.time() - t0
            f.close()
            f = open("pickle.tmp", "rb")
            t0 = time.time()
            b = pickle.load(f)
            tl = time.time() - t0
            f.close()
          if (verbose):
            print array_type.__name__, n, pickler_name, "%.2f %.2f" % (td, tl)
          sys.stdout.flush()
  finally:
    try: os.unlink("pickle.tmp")
    except OSError: pass

def exercise_py_object():
  a = flex.py_object(flex.grid(2,3), value=3)
  assert a.accessor().focus() == (2,3)
  assert a.data() == [3,3,3,3,3,3]
  a = flex.py_object(flex.grid(2,3), value_factory=list)
  assert a.data() == [[],[],[],[],[],[]]
  a = flex.py_object(flex.grid(2,3), values=range(6))
  assert a.data() == range(6)
  assert a[(1,2)] == 5
  a[(1,2)] = -5
  assert a[(1,2)] == -5
  a[(0,0)] = -10
  a[(0,1)] = -1
  a[(0,2)] = -2
  a[(1,0)] = -3
  a[(1,1)] = -4
  assert a.data() == [-10,-1,-2,-3,-4,-5]

def run(iterations):
  i = 0
  while (iterations == 0 or i < iterations):
    exercise_flex_grid()
    exercise_flex_constructors()
    exercise_misc()
    exercise_1d_slicing()
    exercise_push_back_etc()
    exercise_setitem()
    exercise_select()
    exercise_from_stl_vector()
    exercise_operators()
    exercise_bool_inplace_operators()
    exercise_arith_inplace_operators()
    exercise_functions()
    exercise_complex_functions()
    exercise_random()
    exercise_sort()
    exercise_flex_vec3_double()
    exercise_histogram()
    exercise_linear_regression()
    exercise_linear_correlation()
    exercise_mean_and_variance()
    exercise_linear_interpolation()
    exercise_loops()
    exercise_extract_attributes()
    exercise_exceptions()
    exercise_matrix()
    exercise_matrix_cholesky_gill_murray_wright()
    exercise_matrix_move()
    exercise_matrix_inversion_in_place()
    exercise_pickle_single_buffered()
    exercise_pickle_double_buffered()
    pickle_large_arrays(max_exp=2, verbose=0)
    exercise_py_object()
    i += 1

if (__name__ == "__main__"):
  Flags = command_line.parse_options(sys.argv[1:], (
    "large",
  ))
  n = 1
  if (len(sys.argv) > 1 + Flags.n):
    n = int(Flags.regular_args[0])
  if (Flags.large):
    pickle_large_arrays(max_exp=n, verbose=1)
  else:
    run(n)
  print "OK"
