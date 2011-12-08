import scitbx.math
from scitbx.array_family import flex
import math, time, random
from libtbx.test_utils import approx_equal

random.seed(0)
flex.set_random_seed(0)

def residual(a,b,c,d,x):
  return a*x**3+b*x**2+c*x+d

def abcd(x):
  x1,x2,x3 = x
  a = 1
  b = -x1-x2-x3
  c = x1*x3+x2*x3+x1*x2
  d = -x1*x2*x3
  return a,b,c,d

def sol(x, a,b,c,d, check_roots, eps=1.e-6):
  answer = list(x)
  answer.sort()
  for x in answer:
    assert approx_equal(residual(a=a,b=b,c=c,d=d,x=x), 0.0, eps)
  r = scitbx.math.cubic_equation_real(a=a, b=b, c=c, d=d)
  solution = list(r.x)
  solution.sort()
  for ri in r.residual():
    assert approx_equal(ri, 0.0, eps)
  if(check_roots):
    for i in xrange(3):
      assert approx_equal(answer[i], solution[i], eps)

def exercise1():
  v = 1000
  for i in xrange(10000):
    for scale in [0.0000001, 0.0001, 1.]:
      ri1 = random.randint(-v,v)*scale
      ri2 = random.randint(-v,v)*scale
      ri3 = random.randint(-v,v)*scale
      ri = flex.double([ri1,ri2,ri3])
      x = flex.double(flex.random_double_point_on_sphere())*ri
      a,b,c,d = abcd(x=x)
      #
      sol(x=x, a=a,b=b,c=c,d=d, check_roots=True)
      #
      sol(x=(x[0],x[0],x[0]), a=a,b=b,c=c,d=d, check_roots=False)
      sol(x=(x[1],x[1],x[1]), a=a,b=b,c=c,d=d, check_roots=False)
      sol(x=(x[2],x[2],x[2]), a=a,b=b,c=c,d=d, check_roots=False)
      #
      sol(x=(x[0],x[0],x[2]), a=a,b=b,c=c,d=d, check_roots=False)
      sol(x=(x[0],x[2],x[2]), a=a,b=b,c=c,d=d, check_roots=False)
      sol(x=(x[1],x[0],x[1]), a=a,b=b,c=c,d=d, check_roots=False)
      sol(x=(x[1],x[1],x[2]), a=a,b=b,c=c,d=d, check_roots=False)
      sol(x=(x[2],x[0],x[2]), a=a,b=b,c=c,d=d, check_roots=False)

def exercise2():
  a,b,c,d = 3,-10,14,27
  r = scitbx.math.cubic_equation_real(a=a, b=b, c=c, d=d)
  assert approx_equal(r.x[0],-1)
  assert approx_equal(r.x[1], 0) # dummy roots, in reality they are imaginary
  assert approx_equal(r.x[2], 0) # dummy roots, in reality they are imaginary
  assert approx_equal(residual(a=a,b=b,c=c,d=d,x=-1), 0)

def exercise3():
  a,b,c,d = 1,-9,27,-27
  r = scitbx.math.cubic_equation_real(a=a, b=b, c=c, d=d)
  assert approx_equal(r.x[0], 3)
  assert approx_equal(r.x[1], 3)
  assert approx_equal(r.x[2], 3)
  assert approx_equal(residual(a=a,b=b,c=c,d=d,x=3), 0)

def exercise4():
  a,b,c,d = 1,-11,39,-45
  r = scitbx.math.cubic_equation_real(a=a, b=b, c=c, d=d)
  assert approx_equal(r.x[0], 5)
  assert approx_equal(r.x[1], 3)
  assert approx_equal(r.x[2], 3)
  assert approx_equal(residual(a=a,b=b,c=c,d=d,x=3), 0)
  assert approx_equal(residual(a=a,b=b,c=c,d=d,x=5), 0)

def exercise5():
  for i in xrange(100000):
    x1 = random.randint(-10,10)
    x2 = random.randint(-10,10)
    x3 = random.randint(-10,10)
    a,b,c,d = abcd(x = [x1,x2,x3])
    answer = [x1,x2,x3]
    answer.sort()
    answer = flex.double(answer)
    r = scitbx.math.cubic_equation_real(a=a, b=b, c=c, d=d)
    for ri in r.residual():
      assert approx_equal(ri, 0.0)
    solution = list(r.x)
    solution.sort()
    solution = flex.double(solution)
    diff = flex.abs(solution-answer)
    assert approx_equal(diff, [0,0,0])
    for x in [x1,x2,x3, r.x[0],r.x[1],r.x[2]]:
      assert approx_equal(residual(a=a,b=b,c=c,d=d,x=x), 0)

if (__name__ == "__main__"):
  t0 = time.time()
  for i in xrange(1):
    exercise1()
    exercise2()
    exercise3()
    exercise4()
    exercise5()
  print "Time: %6.3f"%(time.time()-t0)
  print "OK"
