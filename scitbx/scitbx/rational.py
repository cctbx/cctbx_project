from scitbx.python_utils import misc
ext = misc.import_ext("scitbx_boost.rational_ext")
misc.import_regular_symbols(globals(), ext.__dict__)
del misc

builtin_int = __builtins__["int"]

def from_string(s):
  flds = [builtin_int(i) for i in s.split("/")]
  assert len(flds) in (1,2)
  if (len(flds) == 1):
    return int(flds[0])
  return int(flds[0], flds[1])

def vector(numerators, denominators):
  if (isinstance(denominators, builtin_int)):
    denominators = [denominators] * len(numerators)
  else:
    assert len(numerators) == len(denominators)
  result = []
  for i in xrange(len(numerators)):
    result.append(int(numerators[i], denominators[i]))
  return result

def lcm_denominators(a):
  l = 1
  for r in a:
    l = lcm(l, r.denominator())
  return l
