
def exercise_real_to_complex_3d (benchmark=True) :
  sizes_1 = [((32,32,32), 16, 0.0000001),
             ((64,64,64), 8, 0.000001),
             ((128,128,128), 8, 0.000001)]#,
# FIXME the 512^3 transform does not yield similar results - why?
#   possibly a FFTW-compatibility issue
#             ((512,512,512), 4, 0.01)]
  sizes_2 = [((36,58,97), 8, 0.000001),
             ((70,120,130), 8, 0.000001),
             ((209,444,320), 4, 0.001),
             ((532,460,485), 4, 0.01)]
  print "real_to_complex_3d"
  _exercise_real_to_complex_3d(sizes_1, benchmark)
  _exercise_real_to_complex_3d(sizes_2, benchmark)

def _exercise_real_to_complex_3d (sizes, benchmark=True) :
  from cctbx import maptbx
  from scitbx.array_family import flex
  from gputbx import cufft
  from scitbx import fftpack
  from libtbx.test_utils import approx_equal
  import time
  for n_real, n_repeats, eps in sizes :
    nx, ny, nz = n_real
    fft = fftpack.real_to_complex_3d((nx,ny,nz))
    mt = flex.mersenne_twister(seed=1)
    g = flex.grid(fft.m_real()).set_focus(fft.n_real())
    map = mt.random_double(size=g.size_1d())
    map.reshape(g)
    sfs = fft.forward(map.deep_copy())
    map2 = fft.backward(sfs)
    sfs_cuda = cufft.real_to_complex_3d_in_place(map)
    map2_cuda = cufft.complex_to_real_3d_in_place(sfs_cuda, n_real)
    maptbx.unpad_in_place(map=map2)
    maptbx.unpad_in_place(map=map2_cuda)
    map2_values = map2.as_1d()
    map2_cuda_values = map2_cuda.as_1d()
    mmm = map2_values.min_max_mean()
    mmm_cuda = map2_cuda_values.min_max_mean()
    assert (map2.size() == map2_cuda.size())
    #assert approx_equal(mmm.min, mmm_cuda.min, eps=0.00001)
    #assert approx_equal(mmm.max, mmm_cuda.max, eps=0.00001)
    assert approx_equal(mmm.mean, mmm_cuda.mean, eps=eps)
    if (benchmark) :
      map_bak = map.deep_copy()
      class c2r_cufft (object) :
        def __init__ (self, n_real) :
          self.n_real = n_real
        def __call__ (self, data) :
          return cufft.complex_to_real_3d_in_place(data, self.n_real)
      cufft_c2r_fn = c2r_cufft(n_real)
      r2c = [fft.forward, cufft.real_to_complex_3d_in_place]
      c2r = [fft.backward, cufft_c2r_fn]
      modules = ["fftpack:", "cufft:  "]
      last_real = [None, None]
      print "  dimensions:", n_real
      print "  repeats:", n_repeats
      k = 0
      for (r2c_fn, c2r_fn, name) in zip(r2c, c2r, modules) :
        t_forward = 0
        t_backward = 0
        map = map_bak.deep_copy()
        for i in range(n_repeats) :
          t1 = time.time()
          sfs = r2c_fn(map)
          t2 = time.time()
          map2 = c2r_fn(sfs)
          t3 = time.time()
          t_forward += t2 - t1
          t_backward += t3 - t2
          if (i == n_repeats - 1) :
            last_real[k] = map2.deep_copy()
        k += 1
        print "    %s %7.3fs (forward)  %7.3fs (backward)" % (name,
          t_forward / n_repeats, t_backward / n_repeats)
      last_fftpack,last_cufft = last_real
      maptbx.unpad_in_place(map=last_fftpack)
      maptbx.unpad_in_place(map=last_cufft)
      mmm = last_fftpack.as_1d().min_max_mean()
      mmm_cuda = last_cufft.as_1d().min_max_mean()
      #assert approx_equal(mmm.min, mmm_cuda.min, eps=0.00001)
      #assert approx_equal(mmm.max, mmm_cuda.max, eps=0.00001)
      assert approx_equal(mmm.mean, mmm_cuda.mean, eps=eps)
      print ""

def exercise_complex_to_complex_3d () :
  from scitbx.array_family import flex
  from gputbx import cufft
  from scitbx import fftpack
  import time
  import sys
  print ""
  print "complex_to_complex_3d"
  for n_complex,n_repeats in [((100,80,90),4), ((200,160,180),4)]:
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
    # XXX extra CuFFT to initialize device
    d = d0.deep_copy()
    cufft.complex_to_complex_3d_in_place(data=d, direction=-1)
    cufft.complex_to_complex_3d_in_place(data=d, direction=+1)
    # benchmarking run
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      cufft.complex_to_complex_3d_in_place(data=d, direction=-1)
      cufft.complex_to_complex_3d_in_place(data=d, direction=+1)
    print "    cufft:    %6.2f seconds" % ((time.time()-t0-overhead)/n_repeats)
    rw = d / np
    #
    t0 = time.time()
    for i_trial in xrange(n_repeats):
      d = d0.deep_copy()
      fftpack.complex_to_complex_3d(n_complex).forward(d)
      fftpack.complex_to_complex_3d(n_complex).backward(d)
    print "    fftpack:  %6.2f seconds" % ((time.time()-t0-overhead)/n_repeats)
    sys.stdout.flush()
    rp = d / np
    #
    print ""
    assert flex.max(flex.abs(rw-rp)) < 1.e-6

if (__name__ == "__main__") :
  exercise_complex_to_complex_3d()
  exercise_real_to_complex_3d()
  print "OK"
