from __future__ import absolute_import, division, print_function
from libtbx.test_utils import approx_equal
from scitbx.array_family import flex
from cma_es import cma_es
from six.moves import range
from six.moves import zip

# =============================================================================
center = flex.double([8.0,-13.0,0.5])

def obj_fun(x):
  assert(len(x) == 3)
  fun = 0.0
  for i in range(len(x)):
    xi = x[i] - center[i]
    fun += xi*xi
  return fun

# from example1.c
def fitfun(x,N):
  isum = 1e4*x[0]*x[0] + 1e-4*x[1]*x[1]
  for i in range(2,N):
    isum += x[i]*x[i]
  return isum

def test_cma_es():
  N = 3
  x = flex.double(N,100.0)
  sd = flex.double(N,5.0)
  m = cma_es(N,x,sd)

  while (not m.converged()):
    # sample population
    p = m.sample_population()
    pop_size = p.accessor().all()[0]

    # update objective function
    v = flex.double(pop_size)
    for i in range(pop_size):
      v[i] = obj_fun(p[(i*N):(i*N + N)])
    m.update_distribution(v)

  x_final = m.get_result()
  assert(approx_equal(x_final,center,eps=1e-6))

def test_cma_es_rosebrock_n(M=10):

  def funct(x,y):
    result = 0
    for xx,yy in zip(x,y):
      result+=100.0*((yy-xx*xx)**2.0) + (1-xx)**2.0
    return result

  N=M*2
  x  = flex.double(N,10.0)
  sd = flex.double(N,3.0)
  m = cma_es(N,x,sd)

  while ( not m.converged() ):
    # sample population
    p = m.sample_population()
    pop_size = p.accessor().all()[0]

    # update objective function
    v = flex.double(pop_size)
    for ii in range(pop_size):
      vector = p[(ii*N):(ii*N + N)]
      x = vector[0:M]
      y = vector[M:]
      v[ii] = funct(x,y)
    m.update_distribution(v)
    print(list(m.get_result()))
    print(flex.min(v))
    print()

  x_final = m.get_result()
  print(list(x_final))

def test_cma_es_lambda():
  N = 3
  x = flex.double(N,100.0)
  sd = flex.double(N,5.0)
  l = 10
  m = cma_es(N,x,sd,l)

  while (not m.converged()):
    # sample population
    p = m.sample_population()
    pop_size = p.accessor().all()[0]

    # update objective function
    v = flex.double(pop_size)
    for i in range(pop_size):
      v[i] = obj_fun(p[(i*N):(i*N + N)])
    m.update_distribution(v)

  x_final = m.get_result()
  assert(approx_equal(x_final,center,eps=1e-6))

def test_cma_es_file():
  import libtbx.load_env
  m = cma_es(libtbx.env.dist_path("cma_es") + "/cma/initials.par")

  while (not m.converged()):
    # sample population and get problem size
    p = m.sample_population()
    pop_size = p.accessor().all()[0]
    N = p.accessor().all()[1]

    # update objective function
    v = flex.double(pop_size)
    for i in range(pop_size):
      v[i] = fitfun(p[(i*N):(i*N + N)],N)
    m.update_distribution(v)

  x_final = m.get_result()
  assert(approx_equal(x_final,flex.double(len(x_final),0.0),eps=1e-5))

# =============================================================================
if (__name__ == '__main__'):
  test_cma_es()
  test_cma_es_lambda()
  test_cma_es_file()
  print('Ok')
