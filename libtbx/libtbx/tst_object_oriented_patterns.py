from libtbx import object_oriented_patterns as oop

def exercise_injector():
  class a(object):
    def __init__(self, i): self.i = i
    def get(self): return self.i
    def set(self, i): self.i = i

  class a_extension(oop.injector, a):
    def get_square(self): return self.get() * self.get()

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
    assert str(err) == "class d has already a method get_square"

def run():
  exercise_injector()
  print 'OK'

if __name__ == '__main__':
  run()
