def exercise():
  import boost.python
  csr = boost.python.ext.string_representation
  from libtbx.str_utils import py_string_representation as psr
  for sr in [csr, psr]:
    assert sr("a", '"', "'") == '"a"'
    assert sr("b", "'", '"') == "'b'"
  def check(s):
    c = csr(s, '"', "'")
    p = psr(s, '"', "'")
    assert c == p
    r = eval(c)
    assert r == s
  iset = range(130) + range(250,256)
  for i in iset:
    s = chr(i)
    check(s)
    for j in iset:
      t = s + chr(j)
      check(t)

def run(args):
  assert args in [[], ["--forever"]]
  while True:
    exercise()
    if (len(args) == 0):
      break
  print "OK"

if (__name__ == "__main__"):
  import sys
  run(args=sys.argv[1:])
