import boost.python
ext = boost.python.import_ext("scitbx_optional_copy_ext")
from scitbx.array_family import shared

def exercise():

  a = shared.unsigned((1, 2, 4))
  t = ext.test(a)
  assert tuple(t.a) == tuple(a)
  assert t.b is None

def run():
  exercise()
  print 'OK'

if __name__ == '__main__':
  run()
