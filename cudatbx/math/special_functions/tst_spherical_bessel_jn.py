import time

from cudatbx import number_of_gpus, reset_gpu
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex

# =============================================================================
def spherical_bessel_jn_test(write_output = False):
  from scitbx.math import spherical_bessel_array
  from cudatbx.math.special_functions import cuda_spherical_bessel_jn

  gpu_id = 0
  z_size = 10000
  z_max = 200.0
  order = 50
  z = flex.double(z_size)
  for i in xrange(z_size):
    z[i] = z_max * (i+1)/z_size

  dt = [0.0,0.0]

  # GPU
  t0 = time.time()
  jn_gpu = cuda_spherical_bessel_jn(order,z,gpu_id)
  t1 = time.time()
  dt[0] = t1 - t0
  if write_output:
    f = open('jn_gpu.dat','w')
    for i in xrange(order+1):
      for j in xrange(z_size):
        f.write('%f %f\n'%(z[j],jn_gpu[i*z_size + j]))
      f.write('&\n')
    f.close()

  # CPU
  jn_cpu = [ None for i in xrange(order+1) ]
  t0 = time.time()
  for n in xrange(order+1):
    jn_cpu[n] = spherical_bessel_array(n,z)
  t1 = time.time()
  dt[1] = t1 - t0
  if write_output:
    f = open('jn_cpu.dat','w')
    for i in xrange(order+1):
      for j in xrange(z_size):
        f.write('%f %f\n'%(z[j],jn_cpu[i][j]))
      f.write('&\n')
    f.close()

  # difference
  d_jn = [ None for i in xrange(order+1) ]
  for n in xrange(order+1):
    d_jn[n] = jn_cpu[n] - jn_gpu[n*z_size:n*z_size + z_size]
    for i in xrange(z_size):
      assert( approx_equal(d_jn[n][i]*d_jn[n][i],0.0,eps=1.0e-6) )
  if write_output:
    f = open('d_jn.dat','w')
    for i in xrange(order+1):
      for j in xrange(z_size):
        f.write('%f %f\n'%(z[j],d_jn[i][j]))
      f.write('&\n')
    f.close()

  return dt

# =============================================================================
if (__name__ == '__main__'):
  import libtbx.load_env
  if (libtbx.env.build_options.enable_cuda):
    t = spherical_bessel_jn_test()

    n_gpus = number_of_gpus()
    for i in xrange(n_gpus):
      reset_gpu(i)

  print 'Ok'
