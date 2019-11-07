from __future__ import absolute_import, division, print_function
import time, random
from cctbx import maptbx
from scitbx.array_family import flex
from libtbx.test_utils import approx_equal
from libtbx import group_args
from six.moves import range

def func(f,x,y,z):
  return eval(f)

def get_exact(f, gx,gy,gz, x,y,z):
  return group_args(
    t  = func(f=f,  x=x, y=y, z=z),
    gx = func(f=gx, x=x, y=y, z=z),
    gy = func(f=gy, x=x, y=y, z=z),
    gz = func(f=gz, x=x, y=y, z=z))

def re(x, xe):
  return abs(x-xe)/abs(xe)*100.

def get_interpolated(site_frac, O):
  # linear
  l0=O.map_data.eight_point_interpolation(site_frac)
  l =O.map_data.eight_point_interpolation_with_gradients(site_frac, O.step)
  assert approx_equal(l[0],l0)
  # quadratic
  q=O.map_data.quadratic_interpolation_with_gradients(site_frac, O.step)
  # tricubic
  t0 = O.map_data.tricubic_interpolation(site_frac)
  t  = O.map_data.tricubic_interpolation_with_gradients(site_frac, O.step)
  # Holds only of "Other option" is used (see maptbx/interpolation.h).
  # assert approx_equal(t[0],t0)
  return group_args(l = l, q = q, t = t)

def set_map(f, n, box=10.):
  n_real = (n,n,n)
  step = box/n
  map_data  = flex.double(flex.grid([n,n,n]))
  for i in range(n):
    x  = i*step
    for j in range(n):
      y  = j*step
      for k in range(n):
        z = k*step
        map_data[i,j,k] = func(f, x=x,y=y,z=z)
  return group_args(
    box=box, step=[step,]*3, map_data=map_data)

def exercise(f,gx,gy,gz,function_type, eps = 1.e-6):
  """
  Exercise linear, quadratic and tricubic interpolation. Thorough.
  Behaviors at corners and boundares are not exercised.
  """
  for n in [25,50,100]:
    O = set_map(f=f, n=n)
    n_trials = 100
    tl,  tq,  tt  = 0,0,0
    gxl, gxq, gxt = 0,0,0
    gyl, gyq, gyt = 0,0,0
    gzl, gzq, gzt = 0,0,0
    for trial in range(n_trials):
      sx = random.randrange(1,30)/10.*random.choice([-1,0,1])
      sy = random.randrange(1,30)/10.*random.choice([-1,0,1])
      sz = random.randrange(1,30)/10.*random.choice([-1,0,1])
      x,y,z = [5+sx, 5+sy, 5+sz]
      e = get_exact(f=f, gx=gx,gy=gy,gz=gz, x=x,y=y,z=z)
      if(abs(e.t)<eps or abs(e.gx)<eps or abs(e.gy)<eps or abs(e.gz)<eps):
        continue
      site_frac = [x/O.box,y/O.box,z/O.box]
      i = get_interpolated(site_frac=site_frac, O=O)
      tl += re(i.l[0], e.t)
      tq += re(i.q[0], e.t)
      tt += re(i.t[0], e.t)
      #
      gxl += re(i.l[1], e.gx)
      gyl += re(i.l[2], e.gy)
      gzl += re(i.l[3], e.gz)
      #
      gxq += re(i.q[1], e.gx)
      gyq += re(i.q[2], e.gy)
      gzq += re(i.q[3], e.gz)
      #
      gxt += re(i.t[1], e.gx)
      gyt += re(i.t[2], e.gy)
      gzt += re(i.t[3], e.gz)
    r = []
    for it in [tl,gxl,gyl,gzl, tq,gxq,gyq,gzq, tt,gxt,gyt,gzt]:
      r.append(it/n_trials)
    tl,gxl,gyl,gzl, tq,gxq,gyq,gzq, tt,gxt,gyt,gzt = r
    #
    s = "(%10.6f %10.6f %10.6f %10.6f)"
    fmt = " ".join(["%5.3f", s, s, s])
    print(fmt%(O.step[0], tl,gxl,gyl,gzl, tq,gxq,gyq,gzq, tt,gxt,gyt,gzt))
    #
    cntr=0
    if(function_type=="L"):
      cntr+=1
      for it in [tl,gxl,gyl,gzl, tq,gxq,gyq,gzq, tt,gxt,gyt,gzt]:
        assert approx_equal(it, 0.0, 1.e-4)
    if(function_type=="Q"):
      cntr+=1
      for it in [tq,gxq,gyq,gzq, tt,gxt,gyt,gzt]:
        assert approx_equal(it, 0.0, 1.e-4)
    if(function_type=="C"):
      cntr+=1
      for it in [tt,gxt,gyt,gzt]:
        assert approx_equal(it, 0.0, 1.e-4)
    assert cntr==1

if (__name__ == "__main__"):
  t0 = time.time()
  fgs = [
  ("L", "3*(x-5)**1-2*(y-5)**1+6*(z-5)**1", "3*(x-5)**0","-2*(y-5)**0", "6*(z-5)**0"),
  ("Q", "3*(x-5)**2-2*(y-5)**2+6*(z-5)**2", "6*(x-5)**1","-4*(y-5)**1","12*(z-5)**1"),
  ("C", "3*(x-5)**3-2*(y-5)**3+6*(z-5)**3", "9*(x-5)**2","-6*(y-5)**2","18*(z-5)**2"),
  ]
  for it in fgs:
    function_type, f,gx,gy,gz = it
    print("- using: f=%s, gx=%s, gy=%s, gz=%s"%(f, gx,gy,gz))
    print(" Step                   Linear                                 "\
          "        Quadratic                              Tricubic")
    exercise(f = f, gx = gx, gy = gy, gz = gz, function_type=function_type)
    print()
  print("Time: %6.2f"%(time.time()-t0))
