from scitbx import math
from scitbx.array_family import flex
from scitbx.array_family import shared
from stdlib import math as smath


def tst_nlm():
  #
  nlm_array = math.nlm_array(10)
  nlm_array.set_coef(10,2,2, 3.0)
  a = nlm_array.get_coef(10,2,2)

  assert ( abs(a-3.0) ) <= 1e-5
  nlm_ind = shared.tiny_int_3( [(10,2,2),(8,2,2)]  )
  nlm_coef = flex.double( [15.0,15.0] )
  assert nlm_array.load_coefs( nlm_ind , nlm_coef)
  a = nlm_array.get_coef(10, 2, 2)
  b = nlm_array.get_coef(8, 2, 2)
  assert ( abs(a-15.0) ) <= 1e-5
  assert ( abs(b-15.0) ) <= 1e-5
  nlm = nlm_array.nlm()
  cnlm = nlm_array.coefs()

  sel = nlm_array.select_on_nl(2,2)

  assert len(sel)==3
  assert 4 in sel
  assert 5 in sel
  assert 6 in sel



def tst_log_factorial():
  # Log factorial generator
  N=89
  lfg =  math.log_factorial_generator(N)
  ff = 1.0
  for ii in range(1,N):
    tmp = lfg.log_fact(ii)
    tmp2 = lfg.fact(ii)
    ff = ff*ii
    #print ii, ff, smath.exp(tmp), tmp2
    assert 1e10*abs(tmp2-ff)/ff<1


def tst_zernike_radial():
  N=50
  M=25
  lfg =  math.log_factorial_generator(N)
  for n in range(M):
    for l in range(n):
      if (n-l)%2==0:
        rzfa = math.zernike_radial(n,l, lfg)
        r = flex.double( flex.double(range(10000))/9999.0)
        a = rzfa.f( r )
        tmp = a*a*r*r
        tmp = flex.sum( tmp )/10000.0
        assert abs(tmp-1)<2e-2
        for nn in range(M):
          for ll in range(nn):
            if (nn-ll)%2==0:
              if (nn!=n):
                if ll==l:
                  rzfb = math.zernike_radial(nn,ll, lfg)
                  b = rzfb.f( r )
                  tmp = a*b*r*r
                  tmp = flex.sum( tmp )/10000.0
                  assert (tmp<1e-2)



def tst_zernike_grid():
  M=10
  N=5
  zg = math.zernike_grid(M,N)
  xyz = zg.xyz()


def tst_zernike_polynome():
  N =50
  n=8
  l=6
  m=6
  lfg =  math.log_factorial_generator(N)
  rnl = math.zernike_radial(n,l,lfg)
  zpl = math.zernike_polynome(n,l,m,rnl)



if __name__ == "__main__":
  tst_nlm()
  tst_log_factorial()
  tst_zernike_radial()
  tst_zernike_grid()
  tst_zernike_polynome()
  print "OK"
