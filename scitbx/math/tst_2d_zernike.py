from scitbx import math
from scitbx.array_family import flex
from scitbx.array_family import shared
from stdlib import math as smath

def tst_2d_zernike_mom(n,l):
  N=20
  nmax = 6
  lfg =  math.log_factorial_generator(N)
  rzfa = math.zernike_2d_radial(n,l,lfg)
  rap = math.zernike_2d_polynome(n,l,rzfa)

  image = flex.vec3_double()

  count = 0
  for x in range(-N, N+1):
    for y in range(-N, N+1):
      rr = smath.sqrt(x*x+y*y)/N
      tt = smath.atan2(y,x)
      value = rap.f(rr,tt)
      value = abs( value )
      if rr>1.0:
        value=0.0
      else:
        count = count + 1
      image.append([x+N,y+N,value])
      print x+N,y+N, value

  grid_2d = math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  zernike_2d_mom = math.two_d_zernike_moments( grid_2d, nmax )

  moments = zernike_2d_mom.moments()

  coefs = moments
  nls = math.nl_array( nmax ).nl()

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
      print x,y,value
      i=i+1

  for nl, c in zip( nls, coefs):
    print nl, c



if __name__ == "__main__":
  tst_2d_zernike_mom(3,1)
  print "OK"
