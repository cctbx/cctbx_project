import fftbx

def show_cseq(vc):
  for i in xrange(len(vc)):
    print "(%.6g,%.6g)" % (vc[i].real, vc[i].imag)
  print

def show_rseq(vr, N):
  for i in xrange(N):
    print "%.6g" % (vr[i],)
  print

def show_rseq_3d(map, N):
  for i in xrange(N[0]):
    for j in xrange(N[1]):
      for k in xrange(N[2]):
        print "%.6g" % (map[(i, j, k)],)
  print

def test_complex_to_complex():
  fft = fftbx.complex_to_complex(10)
  vc = fftbx.vector_of_complex(fft.N())
  for i in xrange(fft.N()):
    vc[i] = complex(2.*i, 2.*i+1.)
  fft.forward(vc)
  show_cseq(vc)
  fft.backward(vc)
  show_cseq(vc)

def test_complex_to_complex_3d():
  fft = fftbx.complex_to_complex_3d((3,4,5))
  vc = fftbx.vector_of_complex()
  map = fftbx.vc3d_accessor(fft.N(), vc, 1)
  map[(0,0,0)] = 123+4j
  assert map[(0,0,0)] == 123+4j
  for i in xrange(vc.size()):
    vc[i] = complex(2.*i, 2.*i+1.)
  fft.forward(map)
  show_cseq(vc)
  fft.backward(map)
  show_cseq(vc)

def test_real_to_complex():
  fft = fftbx.real_to_complex(10)
  vd = fftbx.vector_of_double(fft.Mreal())
  for i in xrange(fft.Nreal()):
    vd[i] = 1.*i
  fft.forward(vd)
  show_rseq(vd, fft.Mreal())
  fft.backward(vd)
  show_rseq(vd, fft.Nreal())

def test_real_to_complex_3d():
  fft = fftbx.real_to_complex_3d((3,4,5))
  vd = fftbx.vector_of_double()
  map = fftbx.vd3d_accessor(fft.Mreal(), vd, 1)
  map[(0,0,0)] = 123
  assert map[(0,0,0)] == 123
  for i in xrange(vd.size()):
    vd[i] = 1.*i
  fft.forward(map)
  show_rseq_3d(map, fft.Mreal())
  fft.backward(map)
  show_rseq_3d(map, fft.Nreal())

if (__name__ == "__main__"):
  test_complex_to_complex()
  test_real_to_complex()
  test_complex_to_complex_3d()
  test_real_to_complex_3d()
