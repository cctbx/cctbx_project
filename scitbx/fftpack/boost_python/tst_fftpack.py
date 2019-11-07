from __future__ import absolute_import, division, print_function
from scitbx import fftpack
from scitbx.array_family import flex
import omptbx # initializes OpenMP environment
from libtbx.test_utils import approx_equal
import libtbx.utils
import random
import math
import sys
from six.moves import range

def fmtfloat(f):
  s = "%.1f" % (f,)
  if (s == "-0.0"): return s[1:]
  return s

def show_cseq(vc):
  for i in range(len(vc)):
    print("(%s,%s)" % (fmtfloat(vc[i].real), fmtfloat(vc[i].imag)))
  print()

def show_rseq(vr, n):
  for i in range(n):
    print(fmtfloat(vr[i]))
  print()

def show_rseq_3d(a, m, n):
  for i in range(n[0]):
    for j in range(n[1]):
      for k in range(n[2]):
        print(fmtfloat(a[(i * m[1] + j) * m[2] + k]))
  print()

def assert_complex_eq_real(vc, vd):
  for i in range(vc.size()):
    assert vc[i].real == vd[2*i]
    assert vc[i].imag == vd[2*i+1]

def test_complex_to_complex(verbose):
  fft = fftpack.complex_to_complex(5)
  vc = flex.complex_double(fft.n())
  vd = flex.double(fft.n() * 2)
  for i in range(fft.n()):
    vc[i] = complex(2.*i, 2.*i+1.)
    vd[2*i] = 2.*i
    vd[2*i+1] = 2.*i+1.
  vct = fft.forward(vc)
  vdt = fft.forward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,)
    assert t.all()[0] == fft.n()
    assert t.focus()[0] == fft.n()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  vct = fft.backward(vc)
  vdt = fft.backward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,)
    assert t.all()[0] == fft.n()
    assert t.focus()[0] == fft.n()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_complex_to_complex_2d(verbose):
  fft = fftpack.complex_to_complex_2d((6,10))
  n = fft.n()
  vc = flex.complex_double(flex.grid(n))
  vd = flex.double(flex.grid((n[0], 2*n[1])))
  for i in range(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  for i in range(vd.size()):
    vd[i] = i
  vct = fft.forward(vc)
  vdt = fft.forward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,0)
    assert t.all() == fft.n()
    assert t.focus() == fft.n()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  vct = fft.backward(vc)
  vdt = fft.backward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,0)
    assert t.all() == fft.n()
    assert t.focus() == fft.n()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  s = vd.size() // 2
  for i in range(vd.size()):
    assert approx_equal(vd[i], s*i)

def test_complex_to_complex_3d(verbose):
  fft = fftpack.complex_to_complex_3d((3,4,5))
  n = fft.n()
  vc = flex.complex_double(flex.grid(n))
  vd = flex.double(flex.grid((n[0], n[1], 2 * n[2])))
  for i in range(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  for i in range(vd.size()):
    vd[i] = i
  vct = fft.forward(vc)
  vdt = fft.forward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,0,0)
    assert t.all() == fft.n()
    assert t.focus() == fft.n()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  vct = fft.backward(vc)
  vdt = fft.backward(vd)
  for t in (vct, vdt):
    assert t.origin() == (0,0,0)
    assert t.all() == fft.n()
    assert t.focus() == fft.n()
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_real_to_complex(verbose):
  fft = fftpack.real_to_complex(6)
  vd = flex.double(fft.m_real())
  vc = flex.complex_double(fft.n_complex())
  for i in range(fft.n_real()):
    vd[i] = 1.*i
  for i in range(fft.n_complex()):
    vc[i] = complex(vd[2*i], vd[2*i+1])
  vdt = fft.forward(vd)
  vct = fft.forward(vc)
  for t in (vdt, vct):
    assert t.origin() == (0,)
    assert t.all()[0] == fft.n_complex()
    assert t.focus()[0] == fft.n_complex()
  if (verbose): show_rseq(vd, fft.m_real())
  assert_complex_eq_real(vc, vd)
  vdt = fft.backward(vd)
  vct = fft.backward(vc)
  for t in (vdt, vct):
    assert t.origin() == (0,)
    assert t.all()[0] == fft.m_real()
    assert t.focus()[0] == fft.n_real()
  if (verbose): show_rseq(vd, fft.n_real())
  assert_complex_eq_real(vc, vd)

  """ The backward fft computes
  y_k = sum_{n=0}^15 x_n e^{i 2pi n/16 k}
      = sum_{n=0}^7 2 Re( x_n e^{i 2pi n/16 k} ) + x_8 cos{i k pi}
  because of the assumed Hermitian condition x_{16-i} = x_i^* (=> x_8 real).
  Only x_0, ..., x_8 need to be stored.
  """
  fft = fftpack.real_to_complex(16)
  assert fft.n_complex() == 9
  x = flex.complex_double(fft.n_complex(), 0)
  x[4] = 1 + 1.j
  y = fft.backward(x)
  for i in range(0, fft.n_real(), 4):
    assert tuple(y[i:i+4]) == (2, -2, -2, 2)

  x = flex.complex_double(fft.n_complex(), 0)
  x[2] = 1 + 1.j
  y = fft.backward(x)
  for i in range(0, fft.n_real(), 8):
    assert tuple(y[i:i+3]) == (2, 0, -2)
    assert approx_equal(y[i+3], -2*math.sqrt(2))
    assert tuple(y[i+4:i+7]) == (-2, 0, 2)
    assert approx_equal(y[i+7], 2*math.sqrt(2))

  x = flex.complex_double(fft.n_complex(), 0)
  x[8] = 1.
  y = fft.backward(x)
  for i in range(0, fft.n_real(), 2):
    assert tuple(y[i:i+2]) == (1, -1)

def exercise_real_to_complex_padding_area():
  mt = flex.mersenne_twister(seed=1)
  for n_real in range(1,101):
    fft = fftpack.real_to_complex(n_real)
    assert fft.n_real() == n_real
    assert fft.n_complex() == n_real // 2 + 1
    assert fft.m_real() == fft.n_complex() * 2
    m_real = fft.m_real()
    z = mt.random_double(size=m_real)*2-1 # non-zero values in padded area
    c = z.deep_copy()
    fft.forward(c)
    if (n_real % 2 == 0):
      assert approx_equal(c[-1], 0) # imaginary part of last complex value
    r = c.deep_copy()
    fft.backward(r)
    r *= (1/n_real)
    assert approx_equal(r[:n_real], z[:n_real])
    if (n_real % 2 == 0):
      assert approx_equal(r[n_real:], [0,0])
      q = c.deep_copy()
      q[-1] = 123 # random imaginary part (which should be zero)
      fft.backward(q)
      q *= (1/n_real)
      assert approx_equal(q, r) # obtain zeros in padded area anyway

def test_real_to_complex_3d(verbose):
  for nx in [3,4]:
    for ny in [4,5]:
      for nz in [7,8]:
        fft = fftpack.real_to_complex_3d((nx,ny,nz))
        m = fft.m_real()
        vd = flex.double(flex.grid(m))
        vc = flex.complex_double(flex.grid((m[0], m[1], m[2]//2)))
        assert vd.size() == 2 * vc.size()
        for i in range(vd.size()):
          vd[i] = random.random()*2-1
        vd_orig = vd.deep_copy()
        for i in range(vc.size()):
          vc[i] = complex(vd[2*i], vd[2*i+1])
        vdt = fft.forward(vd)
        vct = fft.forward(vc)
        for t in (vdt, vct):
          assert t.origin() == (0,0,0)
          assert t.all() == fft.n_complex()
          assert t.focus() == fft.n_complex()
        if (verbose): show_rseq_3d(vd, fft.m_real(), fft.m_real())
        assert_complex_eq_real(vc, vd)
        vdt = fft.backward(vd)
        vct = fft.backward(vc)
        for t in (vdt, vct):
          assert t.origin() == (0,0,0)
          assert t.all() == fft.m_real()
          assert t.focus() == fft.n_real()
        if (verbose): show_rseq_3d(vd, fft.m_real(), fft.n_real())
        assert_complex_eq_real(vc, vd)
        f = nx*ny*nz
        for ix in range(nx):
          for iy in range(ny):
            for iz in range(nz):
              assert approx_equal(vdt[(ix,iy,iz)], vd_orig[(ix,iy,iz)]*f)

def compare_vectors(n, m, v_in, v_tr):
  for i in range(m):
    x = v_tr[i] / n
    assert abs(x - v_in[i]) < 1.e-6, "%d/%d %.6g %.6g" % (i, n, v_in[i], x)

def test_comprehensive_cc_1d(max_transform_size):
  for n in range(1,max_transform_size+1):
    fft = fftpack.complex_to_complex(n)
    m = n * 2
    v_in = flex.double()
    for i in range(m):
      v_in.append(random.random())
    for f,b in ((fft.forward, fft.backward), (fft.backward, fft.forward)):
      v_tr = v_in.deep_copy()
      f(v_tr)
      b(v_tr)
      compare_vectors(n, m, v_in, v_tr)

def test_comprehensive_rc_1d(max_transform_size):
  for n in range(1,max_transform_size+1):
    fft = fftpack.real_to_complex(n)
    m = fft.m_real()
    v_in = flex.double()
    for i in range(n):
      v_in.append(random.random())
    for i in range(n, m):
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
  show_times = libtbx.utils.show_times()
  from scitbx.python_utils import command_line
  flags = command_line.parse_options(sys.argv[1:], [
    "RandomSeed",
    "Verbose",
  ])
  if (not flags.RandomSeed): random.seed(0)
  assert fftpack.adjust_gridding(13, 5) == 15
  assert fftpack.adjust_gridding(13, 5, 6) == 18
  assert fftpack.adjust_gridding_triple((13,22,34), 5) == (15,24,36)
  assert fftpack.adjust_gridding_triple((13,22,34), 5, (6,10,8)) == (18,30,40)
  f = fftpack.factorization(30, False)
  assert f.n() == 30
  assert tuple(f.factors()) == (2, 3, 5)
  test_complex_to_complex(flags.Verbose)
  test_real_to_complex(flags.Verbose)
  test_complex_to_complex_2d(flags.Verbose)
  test_complex_to_complex_3d(flags.Verbose)
  test_real_to_complex_3d(flags.Verbose)
  exercise_real_to_complex_padding_area()
  max_transform_size = 300
  test_comprehensive_cc_1d(max_transform_size)
  test_comprehensive_rc_1d(max_transform_size)
  show_times()

if (__name__ == "__main__"):
  run()
