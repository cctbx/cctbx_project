""" Functions useful to write assertions, especially preconditions """

def is_numeric(x):
  try:
    x+1.
    return True
  except KeyboardInterrupt: raise
  except:
    return False


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
  print 'OK'
