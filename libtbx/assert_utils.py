""" Functions useful to write assertions, especially preconditions """

def is_numeric(x):
  try: x+1.
  except KeyboardInterrupt: raise
  except: return False
  else: return True

def is_string(s):
  try: s+""
  except KeyboardInterrupt: raise
  except: return False
  else: return True

class shall_raise(object):

  def __init__(self, func, *exceptions):
    self.func = func
    self.exceptions = exceptions

  def __call__(self, *args, **kwds):
    try: self.func(*args, **kwds)
    except KeyboardInterrupt: raise
    except:
      import sys
      return sys.exc_info()[0] in self.exceptions
    else: return False


if __name__ == '__main__':
  assert is_numeric(1)
  assert is_numeric(1.)
  assert not is_numeric("1")
  assert not is_numeric([1])
  assert not is_numeric((1,))
  try:
    from scitbx.array_family import flex
  except ImportError:
    import sys
    print >> sys.stderr, "scitbx library not available:" \
                         "some tests were skipped"
  else:
    assert is_numeric(flex.double((1,2,3)))

  def f(x,y,z):
    z//(x-y)
  assert shall_raise(f, ZeroDivisionError)(1,1,1)
  assert not shall_raise(f, ZeroDivisionError)(x=1, y=0, z=0)

  print 'OK'
