import sys, os, random
from scitbx.array_family import flex
from scitbx import fftpack

def fmtfloat(f):
  s = "%.1f" % (f,)
  if (s == "-0.0"): return s[1:]
  return s

def show_cseq(vc):
  for i in xrange(len(vc)):
    print "(%s,%s)" % (fmtfloat(vc[i].real), fmtfloat(vc[i].imag))
  print

def show_rseq(vr, N):
  for i in xrange(N):
    print fmtfloat(vr[i])
  print

def show_rseq_3d(a, M, N):
  for i in xrange(N[0]):
    for j in xrange(N[1]):
      for k in xrange(N[2]):
        print fmtfloat(a[(i * M[1] + j) * M[2] + k])
  print

def assert_complex_eq_real(vc, vd):
  for i in xrange(vc.size()):
    assert vc[i].real == vd[2*i]
    assert vc[i].imag == vd[2*i+1]

def test_complex_to_complex(verbose):
  fft = fftpack.complex_to_complex(5)
  vc = flex.complex_double(fft.N())
  vd = flex.double(fft.N() * 2)
  for i in xrange(fft.N()):
    vc[i] = complex(2.*i, 2.*i+1.)
    vd[2*i] = 2.*i
    vd[2*i+1] = 2.*i+1.
  vct = fft.forward(vc)
  vdt = fft.forward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,)
    assert t.grid()[0] == fft.N()
    assert t.layout()[0] == fft.N()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  vct = fft.backward(vc)
  vdt = fft.backward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,)
    assert t.grid()[0] == fft.N()
    assert t.layout()[0] == fft.N()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_complex_to_complex_3d(verbose):
  fft = fftpack.complex_to_complex_3d((3,4,5))
  N = fft.N()
  vc = flex.complex_double(flex.grid(N))
  vd = flex.double(flex.grid((N[0], N[1], 2 * N[2])))
  for i in xrange(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  for i in xrange(vd.size()):
    vd[i] = i
  vct = fft.forward(vc)
  vdt = fft.forward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,0,0)
    assert t.grid() == fft.N()
    assert t.layout() == fft.N()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  vct = fft.backward(vc)
  vdt = fft.backward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,0,0)
    assert t.grid() == fft.N()
    assert t.layout() == fft.N()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_real_to_complex(verbose):
  fft = fftpack.real_to_complex(6)
  vd = flex.double(fft.Mreal())
  vc = flex.complex_double(fft.Ncomplex())
  for i in xrange(fft.Nreal()):
    vd[i] = 1.*i
  for i in xrange(fft.Ncomplex()):
    vc[i] = complex(vd[2*i], vd[2*i+1])
  vdt = fft.forward(vd)
  vct = fft.forward(vc)
  for t in (vdt, vct):
    assert t.origin() == (0,)
    assert t.grid()[0] == fft.Ncomplex()
    assert t.layout()[0] == fft.Ncomplex()
  if (verbose): show_rseq(vd, fft.Mreal())
  assert_complex_eq_real(vc, vd)
  vdt = fft.backward(vd)
  vct = fft.backward(vc)
  for t in (vdt, vct):
    assert t.origin() == (0,)
    assert t.grid()[0] == fft.Mreal()
    assert t.layout()[0] == fft.Nreal()
  if (verbose): show_rseq(vd, fft.Nreal())
  assert_complex_eq_real(vc, vd)

def test_real_to_complex_3d(verbose):
  fft = fftpack.real_to_complex_3d((3,4,5))
  M = fft.Mreal()
  vd = flex.double(flex.grid(M))
  vc = flex.complex_double(flex.grid((M[0], M[1], M[2]/2)))
  assert vd.size() == 2 * vc.size()
  for i in xrange(vd.size()):
    vd[i] = i
  for i in xrange(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  vdt = fft.forward(vd)
  vct = fft.forward(vc)
  for t in (vdt, vct):
    assert t.origin() == (0,0,0)
    assert t.grid() == fft.Ncomplex()
    assert t.layout() == fft.Ncomplex()
  if (verbose): show_rseq_3d(vd, fft.Mreal(), fft.Mreal())
  assert_complex_eq_real(vc, vd)
  vdt = fft.backward(vd)
  vct = fft.backward(vc)
  for t in (vdt, vct):
    assert t.origin() == (0,0,0)
    assert t.grid() == fft.Mreal()
    assert t.layout() == fft.Nreal()
  if (verbose): show_rseq_3d(vd, fft.Mreal(), fft.Nreal())
  assert_complex_eq_real(vc, vd)

def compare_vectors(n, m, v_in, v_tr):
  for i in xrange(m):
    x = v_tr[i] / n
    assert abs(x - v_in[i]) < 1.e-6, "%d/%d %.6g %.6g" % (i, n, v_in[i], x)

def test_comprehensive_cc_1d(max_transform_size):
  for n in xrange(1,max_transform_size+1):
    fft = fftpack.complex_to_complex(n)
    m = n * 2
    v_in = flex.double()
    for i in xrange(m):
      v_in.append(random.random())
    for f,b in ((fft.forward, fft.backward), (fft.backward, fft.forward)):
      v_tr = v_in.deep_copy()
      f(v_tr)
      b(v_tr)
      compare_vectors(n, m, v_in, v_tr)

def test_comprehensive_rc_1d(max_transform_size):
  for n in xrange(1,max_transform_size+1):
    fft = fftpack.real_to_complex(n)
    m = fft.Mreal()
    v_in = flex.double()
    for i in xrange(n):
      v_in.append(random.random())
    for i in xrange(n, m):
      v_in.append(999)
    v_tr = v_in.deep_copy()
    fft.forward(v_tr)
    fft.backward(v_tr)
    compare_vectors(n, n, v_in, v_tr)
    v_in[n] = v_in[1]
    v_in[1] = 0
    if (n % 2 == 0): v_in[n+1] = 0
    v_tr = v_in.deep_copy()
    fft.backward(v_tr)
    fft.forward(v_tr)
    compare_vectors(n, m, v_in, v_tr)

def run():
  from scitbx.python_utils import command_line
  Flags = command_line.parse_options(sys.argv[1:], (
    "RandomSeed",
    "verbose",
  ))
  if (not Flags.RandomSeed): random.seed(0)
  assert fftpack.adjust_gridding(13, 5) == 15
  assert fftpack.adjust_gridding(13, 5, 6) == 18
  assert fftpack.adjust_gridding_triple((13,22,34), 5) == (15,24,36)
  assert fftpack.adjust_gridding_triple((13,22,34), 5, (6,10,8)) == (18,30,40)
  f = fftpack.factorization(30, 0)
  assert f.N() == 30
  assert tuple(f.Factors()) == (2, 3, 5)
  test_complex_to_complex(Flags.verbose)
  test_real_to_complex(Flags.verbose)
  test_complex_to_complex_3d(Flags.verbose)
  test_real_to_complex_3d(Flags.verbose)
  max_transform_size = 300
  test_comprehensive_cc_1d(max_transform_size)
  test_comprehensive_rc_1d(max_transform_size)

if (__name__ == "__main__"):
  run()
  t = os.times()
  print "u+s,u,s: %.2f %.2f %.2f" % (t[0] + t[1], t[0], t[1])
