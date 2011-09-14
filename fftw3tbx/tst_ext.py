import fftw3tbx
import omptbx # initializes OpenMP environment
from scitbx import fftpack
from scitbx.array_family import flex
import libtbx.utils
import time
import sys

try: from cctbx import maptbx
except ImportError: maptbx = None

def exercise_complex_to_complex():
  print "complex_to_complex"
  for n in xrange(1,256+1):
    dp = flex.polar(
      flex.random_double(size=n)*2-1,
      flex.random_double(size=n)*2-1)
    dw = dp.deep_copy()
    fft = fftpack.complex_to_complex(n)
    fftw3tbx.complex_to_complex_in_place(data=dw, exp_sign=-1)
    fft.forward(dp)
    assert flex.max(flex.abs(dw-dp)) < 1.e-6
    fftw3tbx.complex_to_complex_in_place(data=dw, exp_sign=+1)
    fft.backward(dp)
    assert flex.max(flex.abs(dw-dp)) < 1.e-6
  for n,n_repeats in [(1200,500), (9600,250)]:
    print "  factors of %d:" % n, list(fftpack.complex_to_complex(n).factors())
    print "  repeats:", n_repeats
    d0 = flex.polar(
      flex.random_double(size=n)*2-1,
      flex.random_double(size=n)*2-1)
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
    overhead = time.time()-t0
    print "    overhead: %.2f seconds" % overhead
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      fftw3tbx.complex_to_complex_in_place(data=d, exp_sign=-1)
      fftw3tbx.complex_to_complex_in_place(data=d, exp_sign=+1)
    print "    fftw:     %.2f seconds" % (time.time()-t0-overhead)
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      fftpack.complex_to_complex(n).forward(d)
      fftpack.complex_to_complex(n).backward(d)
    print "    fftpack:  %.2f seconds" % (time.time()-t0-overhead)
    sys.stdout.flush()

def exercise_complex_to_complex_3d():
  print "complex_to_complex_3d"
  for n_complex,n_repeats in [((100,80,90),2), ((200,160,180),1)]:
    print "  dimensions:", n_complex
    print "  repeats:", n_repeats
    np = n_complex[0]*n_complex[1]*n_complex[2]
    d0 = flex.polar(
      flex.random_double(size=np)*2-1,
      flex.random_double(size=np)*2-1)
    d0.reshape(flex.grid(n_complex))
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
    overhead = time.time()-t0
    print "    overhead: %.2f seconds" % overhead
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      fftw3tbx.complex_to_complex_3d_in_place(data=d, exp_sign=-1)
      fftw3tbx.complex_to_complex_3d_in_place(data=d, exp_sign=+1)
    print "    fftw:     %.2f seconds" % (time.time()-t0-overhead)
    rw = d / np
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      fftpack.complex_to_complex_3d(n_complex).forward(d)
      fftpack.complex_to_complex_3d(n_complex).backward(d)
    print "    fftpack:  %.2f seconds" % (time.time()-t0-overhead)
    sys.stdout.flush()
    rp = d / np
    #
    assert flex.max(flex.abs(rw-rp)) < 1.e-6

def exercise_real_to_complex():
  print "real_to_complex"
  for n in xrange(1,256+1):
    fft = fftpack.real_to_complex(n)
    dp = flex.random_double(size=n)*2-1
    dp.resize(flex.grid(fft.m_real()).set_focus(n))
    dw = dp.deep_copy()
    cw = fftw3tbx.real_to_complex_in_place(dw)
    cp = fft.forward(dp)
    assert flex.max(flex.abs(cw-cp)) < 1.e-6
    rw = fftw3tbx.complex_to_real_in_place(cw, n)
    rp = fft.backward(cp)
    assert flex.max(flex.abs(rw[:n]-rp[:n])) < 1.e-6
  for n,n_repeats in [(2400,500), (19200,250)]:
    fft = fftpack.real_to_complex(n)
    print "  factors of %d:" % n, list(fft.factors())
    print "  repeats:", n_repeats
    d0 = flex.random_double(size=n)*2-1
    d0.resize(flex.grid(fft.m_real()).set_focus(n))
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
    overhead = time.time()-t0
    print "    overhead: %.2f seconds" % overhead
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      c = fftw3tbx.real_to_complex_in_place(data=d)
      fftw3tbx.complex_to_real_in_place(data=c, n=n)
    print "    fftw:     %.2f seconds" % (time.time()-t0-overhead)
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      c = fftpack.real_to_complex(n).forward(d)
      fftpack.real_to_complex(n).backward(c)
    print "    fftpack:  %.2f seconds" % (time.time()-t0-overhead)
    sys.stdout.flush()

def exercise_real_to_complex_3d():
  print "real_to_complex_3d"
  for n_real,n_repeats in [((100,80,90),8),
                           ((200,160,180),2),
                           ((300,240,320),1)]:
    print "  dimensions:", n_real
    print "  repeats:", n_repeats
    fft = fftpack.real_to_complex_3d(n_real)
    m_real = fft.m_real()
    np = n_real[0]*n_real[1]*n_real[2]
    mp = m_real[0]*m_real[1]*m_real[2]
    d0 = flex.random_double(size=mp)*2-1
    d0.reshape(flex.grid(m_real).set_focus(n_real))
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
    overhead = time.time()-t0
    print "    overhead: %.2f seconds" % overhead
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      c = fftw3tbx.real_to_complex_3d_in_place(data=d)
      assert c.all() == fft.n_complex()
      assert c.focus() == fft.n_complex()
      assert c.id() == d.id()
      r = fftw3tbx.complex_to_real_3d_in_place(data=c, n=n_real)
      assert r.all() == fft.m_real()
      assert r.focus() == fft.n_real()
      assert r.id() == d.id()
    print "    fftw:     %.2f seconds" % (time.time()-t0-overhead)
    if (maptbx is not None):
      maptbx.unpad_in_place(map=d)
      rw = d / np
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      c = fftpack.real_to_complex_3d(n_real).forward(d)
      assert c.all() == fft.n_complex()
      assert c.focus() == fft.n_complex()
      assert c.id() == d.id()
      r = fftpack.real_to_complex_3d(n_real).backward(c)
      assert r.all() == fft.m_real()
      assert r.focus() == fft.n_real()
      assert r.id() == d.id()
    print "    fftpack:  %.2f seconds" % (time.time()-t0-overhead)
    sys.stdout.flush()
    if (maptbx is not None):
      maptbx.unpad_in_place(map=d)
      rp = d / np
      #
      assert flex.max(flex.abs(rw-rp)) < 1.e-6

def exercise(args):
  show_times = libtbx.utils.show_times()
  if (fftw3tbx.ext is None):
    print "Skipping fftw3tbx tests: Python extension not available."
  else:
    print "fftw_version:", fftw3tbx.fftw_version
    print "### NOTE: Showing wall-clock times. ###"
    sys.stdout.flush()
    exercise_complex_to_complex()
    exercise_complex_to_complex_3d()
    exercise_real_to_complex()
    exercise_real_to_complex_3d()
  show_times()

if (__name__ == "__main__"):
  exercise(sys.argv[1:])
