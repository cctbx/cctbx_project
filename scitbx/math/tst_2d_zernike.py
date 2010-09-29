from scitbx import math
from scitbx.array_family import flex
from scitbx.array_family import shared
from stdlib import math as smath
import random
import time, os, sys

def read_data(filename):
  file=open(filename, 'r')
  data=flex.vec3_double()
  for line in file:
    keys=line.split()
    if(len(keys)==3):
      x=int(keys[0])
      y=int(keys[1])
      z=float(keys[2])
      data.append([x,y,z])
  file.close()
  return data

def generate_image(n,l, N=100):
  nmax = max(20,n)
  lfg =  math.log_factorial_generator(nmax)
  rzfa = math.zernike_2d_radial(n,l,lfg)
  rap = math.zernike_2d_polynome(n,l,rzfa)

  image = flex.vec3_double()

  original=open('original.dat','w')
  count = 0
  for x in range(-N, N+1):
    for y in range(-N, N+1):
      rr = smath.sqrt(x*x+y*y)/N
      tt = smath.atan2(y,x)
      value = rap.f(rr,tt)
      value = value.real
      if rr>1.0:
        value=0.0
      else:
        count = count + 1
      image.append([x+N,y+N,value])
      print>>original, x+N,y+N, value
  original.close()
  return image

def tst_2d_zernike_mom(n,l, N=100, filename=None):
  nmax = max(20,n)
  rebuilt=open('rebuilt.dat','w')
  tt1=time.time()
  if(filename is not None):
    image=read_data(filename)
  else:
    image=generate_image(n,l)

  NP=int(smath.sqrt( image.size() ))
  N=NP/2
  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  tt2=time.time()
  print "time used: ", tt2-tt1
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )

  moments = zernike_2d_mom.moments()

  tt2=time.time()
  print "time used: ", tt2-tt1

  coefs = flex.real( moments )
  nl_array = math.nl_array( nmax )
  nls = nl_array.nl()
  nl_array.load_coefs( nls, coefs )
  lfg =  math.log_factorial_generator(nmax)

  print nl_array.get_coef(n,l)*2

  for nl, c in zip( nls, moments):
    if(abs(c)<1e-3):
      c=0
    print nl, c

  print
  reconst=flex.complex_double(NP**2, 0)
  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    if(l>0):
      c=c*2
    rzfa = math.zernike_2d_radial(n,l,lfg)
    rap = math.zernike_2d_polynome(n,l,rzfa)
    i=0
    for x in range(0,NP):
      x=x-N
      for y in range(0,NP):
        y=y-N
        rr = smath.sqrt(x*x+y*y)/N
        if rr>1.0:
          value=0.0
        else:
          tt = smath.atan2(y,x)
          value = rap.f(rr,tt)
        reconst[i]=reconst[i]+value*c
        i=i+1

  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      if(value>0):
        print>>rebuilt, x,y,image[i][2],value
      i=i+1


  rebuilt.close()


def tst_2d_poly(n,l):
  nmax=max(n,20)
  np=50
  x,y=0.1,0.9
  r,t=smath.sqrt(x*x+y*y),smath.atan2(y,x)
  lfg =  math.log_factorial_generator(nmax)
  rzfa = math.zernike_2d_radial(n,l,lfg)
  rap = math.zernike_2d_polynome(n,l,rzfa)
  rt_value=rap.f(r,t)
  grid = math.two_d_grid(np, nmax)
  zm2d = math.two_d_zernike_moments(grid, nmax)
  xy_value=zm2d.zernike_poly(n,l,x,y)

  print rt_value, xy_value, abs(rt_value), abs(xy_value)

def tst_2d_zm(n,l):
  nmax=max(n,20)
  np=100
  points=flex.double(range(-np,np+1))/np
  grid = math.two_d_grid(np, nmax)
  zm2d = math.two_d_zernike_moments(grid, nmax)
  image = flex.vec3_double()

  output=file('testmap.dat','w')

  for x in points:
    for y in points:
      r=smath.sqrt(x*x+y*y)
      if(r>1.0):
        value=0.0
      else:
        value=zm2d.zernike_poly(n,l,x,y).real
      image.append([x*np+np,y*np+np, value])

  grid.clean_space( image )
  grid.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid, nmax )

  moments = zernike_2d_mom.moments()

  coefs = flex.real( moments )
  nl_array = math.nl_array( nmax )
  nls = nl_array.nl()

  for nl, c in zip( nls, moments):
    if(abs(c)<1e-3):
      c=0
    print nl, c

  NP=np*2+1
  reconst=zernike_2d_mom.zernike_map(nmax, np)
  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      if(value>0):
        print>>output, x,y,image[i][2],value
      i=i+1


  output.close()



if __name__ == "__main__":
  t1=time.time()
  tst_2d_zm(41,21)
  t2=time.time()
  print"xy:  ", t2-t1
  tst_2d_zernike_mom(41,21)
  t3=time.time()
  print"rt:  ", t3-t2
  exit()

  tst_2d_poly(51,19)
  exit()
  t1 = time.time()
  args = sys.argv[1:]
  filename=None
  if( len(args) > 0 ):
    filename = args[0]
  tst_2d_zernike_mom(30,2, filename=filename)
  exit()

  tst_2d_zernike_mom(17,3)
  tst_2d_zernike_mom(16,2)
  tst_2d_zernike_mom(14,12)
  #tst_2d_random(20)
  t2 = time.time()
  print "time used: ", t2-t1
  print "OK"
