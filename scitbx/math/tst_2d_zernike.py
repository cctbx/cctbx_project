from scitbx import math
from scitbx.array_family import flex
from scitbx.array_family import shared
from stdlib import math as smath
import random
import time

def tst_2d_zernike_mom(n,l):
  N=32
  nmax = n
  lfg =  math.log_factorial_generator(nmax)
  rzfa = math.zernike_2d_radial(n,l,lfg)
  rap = math.zernike_2d_polynome(n,l,rzfa)

  image = flex.vec3_double()

  original=open('original.dat','w')
  rebuilt=open('rebuilt.dat','w')
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


  tt1=time.time()
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
  
  print nl_array.get_coef(n,l)*2

  for nl, c in zip( nls, coefs):
    if(abs(c)<1e-3):
      c=0
    print nl, c

  print
  reconst=flex.complex_double((N*2+1)**2, 0)
  for nl,c in zip( nls, coefs):
    n=nl[0]
    l=nl[1]
    rzfa = math.zernike_2d_radial(n,l,lfg)
    rap = math.zernike_2d_polynome(n,l,rzfa)
    i=0
    for x in range(-N, N+1):
      for y in range(-N, N+1):
        rr = smath.sqrt(x*x+y*y)/N
        if rr>1.0:
          value=0.0
        else:
          tt = smath.atan2(y,x)
          value = rap.f(rr,tt)
        reconst[i]=reconst[i]+value*c
        i=i+1

  i = 0
  for x in range(0,N*2+1):
    for y in range(0,N*2+1):
      value=reconst[i].real
      print>>rebuilt, x,y,value*2
      i=i+1


  original.close()
  rebuilt.close()


if __name__ == "__main__":
<<<<<<< .mine
  t1 = time.time()
  tst_2d_zernike_mom(16,10)
  exit()

  tst_2d_zernike_mom(17,3)
  tst_2d_zernike_mom(16,2)
  tst_2d_zernike_mom(14,12)
  #tst_2d_random(20)
  t2 = time.time()
  print "time used: ", t2-t1
=======
  tst_2d_zernike_mom(3,1)
  print "OK"
>>>>>>> .r11409
