from libtbx import object_oriented_patterns as oop
from libtbx.test_utils import Exception_expected
import libtbx

def exercise_injector():
  class a(object):
    """ doc for a """
    def __init__(self, i): self.i = i
    def get(self): return self.i
    def set(self, i): self.i = i

  class a_extension(oop.injector, a):
    """ doc for extension of a """
    var = 1
    def get_square(self): return self.get() * self.get()

  assert a.var == 1
  o = a(2)
  assert o.get() == 2
  assert o.get_square() == 4

  class b(object):
    def __init__(self, i): self.i = i
    def get(self): return self.i-1

  class c(b):
    def get(self): return self.i + 1

  class c_extension(oop.injector, c):
    def get_square(self): return self.get() * self.get()

  o = c(-3)
  assert o.get() == -2
  assert o.get_square() == 4

  try:
    class d(a): pass
    class d_extension(oop.injector, d):
      def get_square(self): return 0
  except AssertionError, err:
    assert str(err) == "class d already has attribute 'get_square'"

def exercise_memoize():
  diagnostic = []
  def f(x):
    """ Documentation for function f """
    diagnostic.append('+')
    return x+1

  mf = oop.memoize(f)
  assert mf(0) == 1
  assert diagnostic == ['+']
  assert mf(0) == 1
  assert diagnostic == ['+']
  assert mf(1) == 2
  assert diagnostic == ['+']*2
  assert mf(0) == 1
  assert diagnostic == ['+']*2
  assert mf(1) == 2
  assert diagnostic == ['+']*2
  try:
    mf(x=1)
  except TypeError, e:
    pass
  else:
    raise Exception_expected
  assert mf.__doc__ == """ Documentation for function f """


  class foo(object):
    def __init__(self, a):
      self.a = a
    def f(self, x):
      """ Documentation for method f """
      diagnostic.append('+')
      return self.a + x
    f = oop.memoize_method(f)

  diagnostic = []
  o1 = foo(1)
  o2 = foo(2)
  assert o1.f(4) == 5
  assert diagnostic == ['+']
  assert o2.f(5) == 7
  assert diagnostic == ['+']*2
  assert o1.f(2) == 3
  assert diagnostic == ['+']*3
  assert o1.f(4) == 5
  assert diagnostic == ['+']*3
  assert o2.f(5) == 7
  assert diagnostic == ['+']*3
  o3 = o1
  assert o3.f(4) == 5
  assert diagnostic == ['+']*3
  assert o3.f(1) == 2
  assert diagnostic == ['+']*4
  assert o1.f(1) == 2
  assert diagnostic == ['+']*4
  assert o1.f.__doc__ == " Documentation for method f "

  class _foo_extension(oop.injector, foo):
    def g(self, n):
      return 'x'*n
    g = oop.memoize_method(g)

  efoo = foo(0)
  assert efoo.g(2) == "xx"
  assert efoo.g(3) == "xxx"
  assert efoo.g(2) == "xx"
  assert efoo._memoized_g.cached == {(2,):'xx', (3,): 'xxx'}

def exercise_null():
  n = oop.null("a", 1, [])
  assert not n
  assert isinstance(n.f, oop.null)
  assert isinstance(n.g(1, b=2), oop.null)
  h = n.h = 2
  assert h == 2
  assert isinstance(n.h, oop.null)
  assert isinstance(n[23085], oop.null)
  x = n[234] = 2
  assert x == 2
  assert isinstance(n[234], oop.null)

def exercise_easy_property():

  class foo(object):

    def __init__(self):
      self.a = 1

    class bar(libtbx.property):
      """ documentation for property bar """

      def fget(self):
        return self.a

      def fset(self, a):
        self.a = a

      def fdel(self):
        self.a = -1

  x = foo()
  assert x.bar == 1
  x.bar = 2
  assert x.bar == 2
  assert x.a == 2
  x.bar += 4
  assert x.bar == 6
  assert x.a == 6
  assert foo.bar.__doc__.find("documentation for property bar") >= 0
  del x.bar
  assert x.a == -1

def run():
  exercise_injector()
  exercise_memoize()
  exercise_null()
  exercise_easy_property()
  print 'OK'

if __name__ == '__main__':
  run()
