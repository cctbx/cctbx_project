import fftbx

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

def show_rseq_3d(map, N):
  for i in xrange(N[0]):
    for j in xrange(N[1]):
      for k in xrange(N[2]):
        print fmtfloat(map[(i, j, k)])
  print

def assert_complex_eq_real(vc, vd):
  for i in xrange(vc.size()):
    assert vc[i].real == vd[2*i]
    assert vc[i].imag == vd[2*i+1]

def test_complex_to_complex(verbose):
  fft = fftbx.complex_to_complex(5)
  vc = fftbx.vector_of_complex(fft.N())
  vd = fftbx.vector_of_double(fft.N() * 2)
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

def test_complex_to_complex_3d(test_map_cast, verbose):
  fft = fftbx.complex_to_complex_3d((3,4,5))
  vc = fftbx.vector_of_complex()
  mapc = fftbx.vc3d_accessor(fft.N(), vc, 1)
  mapc[(0,0,0)] = 123+4j
  assert mapc[(0,0,0)] == 123+4j
  for i in xrange(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  Ncomplex = fft.N()
  Nreal = (Ncomplex[0], Ncomplex[1], Ncomplex[2] * 2)
  vd = fftbx.vector_of_double()
  mapd = fftbx.vd3d_accessor(Nreal, vd, 1)
  for i in xrange(vd.size()):
    vd[i] = i
  assert vd.size() == 2 * vc.size()
  if (test_map_cast):
    mapc = fftbx.vd3d_accessor(Nreal, vc, 1)
    assert vd.size() == 2 * vc.size()
    mapd = fftbx.vc3d_accessor(fft.N(), vd, 1)
    assert vd.size() == 2 * vc.size()
  fft.forward(mapc)
  fft.forward(mapd)
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)
  fft.backward(mapc)
  fft.backward(mapd)
  if (verbose): show_cseq(vc)
  assert_complex_eq_real(vc, vd)

def test_real_to_complex(verbose):
  fft = fftbx.real_to_complex(6)
  vd = fftbx.vector_of_double(fft.Mreal())
  vc = fftbx.vector_of_complex(fft.Ncomplex())
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

def test_real_to_complex_3d(test_map_cast, verbose):
  fft = fftbx.real_to_complex_3d((3,4,5))
  vd = fftbx.vector_of_double()
  mapd = fftbx.vd3d_accessor(fft.Mreal(), vd, 1)
  mapd[(0,0,0)] = 123
  assert mapd[(0,0,0)] == 123
  for i in xrange(vd.size()):
    vd[i] = i
  vc = fftbx.vector_of_complex()
  mapc = fftbx.vc3d_accessor(fft.Ncomplex(), vc, 1)
  for i in xrange(vc.size()):
    vc[i] = complex(2*i, 2*i+1)
  assert vd.size() == 2 * vc.size()
  if (test_map_cast):
    mapd = fftbx.vc3d_accessor(fft.Ncomplex(), vd, 1)
    assert vd.size() == 2 * vc.size()
    mapc = fftbx.vd3d_accessor(fft.Mreal(), vc, 1)
    assert vd.size() == 2 * vc.size()
  fft.forward(mapd)
  fft.forward(mapc)
  if (verbose): show_rseq_3d(mapd, fft.Mreal())
  assert_complex_eq_real(vc, vd)
  fft.backward(mapd)
  fft.backward(mapc)
  if (verbose): show_rseq_3d(mapd, fft.Nreal())
  assert_complex_eq_real(vc, vd)

if (__name__ == "__main__"):
  import sys
  verbose = not "--quiet" in sys.argv[1:]
  test_complex_to_complex(verbose)
  test_real_to_complex(verbose)
  test_complex_to_complex_3d(0, verbose)
  test_complex_to_complex_3d(1, 0)
  test_real_to_complex_3d(0, verbose)
  test_real_to_complex_3d(1, 0)
