import sys, os
from cctbx_boost.arraytbx import shared
from cctbx_boost import fftbx
from cctbx.development import debug_utils

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
  fft = fftbx.complex_to_complex(5)
  vc = shared.complex_double(fft.N())
  vd = shared.double(fft.N() * 2)
  for i in xrange(fft.N()):
    vc[i] = complex(2.*i, 2.*i+1.)
    vd[2*i] = 2.*i
    vd[2*i+1] = 2.*i+1.
  fft.forward(vc)
  fft.forward(vd)
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  fft.backward(vc)
  fft.backward(vd)
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_complex_to_complex_3d(verbose):
  fft = fftbx.complex_to_complex_3d((3,4,5))
  N = fft.N()
  vc = shared.complex_double(N[0] * N[1] * N[2])
  vd = shared.double(2 * vc.size())
  for i in xrange(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  for i in xrange(vd.size()):
    vd[i] = i
  fft.forward(vc)
  fft.forward(vd)
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  fft.backward(vc)
  fft.backward(vd)
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_real_to_complex(verbose):
  fft = fftbx.real_to_complex(6)
  vd = shared.double(fft.Mreal())
  vc = shared.complex_double(fft.Ncomplex())
  for i in xrange(fft.Nreal()):
    vd[i] = 1.*i
  for i in xrange(fft.Ncomplex()):
    vc[i] = complex(vd[2*i], vd[2*i+1])
  fft.forward(vd)
  fft.forward(vc)
  if (verbose): show_rseq(vd, fft.Mreal())
  assert_complex_eq_real(vc, vd)
  fft.backward(vd)
  fft.backward(vc)
  if (verbose): show_rseq(vd, fft.Nreal())
  assert_complex_eq_real(vc, vd)

def test_real_to_complex_3d(verbose):
  fft = fftbx.real_to_complex_3d((3,4,5))
  M = fft.Mreal()
  vd = shared.double(M[0] * M[1] * M[2])
  vc = shared.complex_double(vd.size() / 2)
  for i in xrange(vd.size()):
    vd[i] = i
  for i in xrange(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  assert vd.size() == 2 * vc.size()
  fft.forward(vd)
  fft.forward(vc)
  if (verbose): show_rseq_3d(vd, fft.Mreal(), fft.Mreal())
  assert_complex_eq_real(vc, vd)
  fft.backward(vd)
  fft.backward(vc)
  if (verbose): show_rseq_3d(vd, fft.Mreal(), fft.Nreal())
  assert_complex_eq_real(vc, vd)

def compare_vectors(n, m, v_in, v_tr):
  for i in xrange(m):
    x = v_tr[i] / n
    assert abs(x - v_in[i]) < 1.e-6, "%d/%d %.6g %.6g" % (i, n, v_in[i], x)

def test_comprehensive_cc_1d(max_transform_size):
  for n in xrange(1,max_transform_size+1):
    fft = fftbx.complex_to_complex(n)
    m = n * 2
    v_in = shared.double()
    for i in xrange(m):
      v_in.append(debug_utils.random.random())
    for f,b in ((fft.forward, fft.backward), (fft.backward, fft.forward)):
      v_tr = v_in.deep_copy()
      f(v_tr)
      b(v_tr)
      compare_vectors(n, m, v_in, v_tr)

def test_comprehensive_rc_1d(max_transform_size):
  for n in xrange(1,max_transform_size+1):
    fft = fftbx.real_to_complex(n)
    m = fft.Mreal()
    v_in = shared.double()
    for i in xrange(n):
      v_in.append(debug_utils.random.random())
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
  Flags = debug_utils.command_line_options(sys.argv[1:], (
    "RandomSeed",
    "verbose",
  ))
  if (not Flags.RandomSeed): debug_utils.set_random_seed(0)
  assert fftbx.adjust_gridding(13, 5) == 15
  assert fftbx.adjust_gridding(13, 5, 6) == 18
  assert fftbx.adjust_gridding_triple((13,22,34), 5) == (15,24,36)
  assert fftbx.adjust_gridding_triple((13,22,34), 5, (6,10,8)) == (18,30,40)
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
  print "u+s,u,s:", t[0] + t[1], t[0], t[1]
