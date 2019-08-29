from __future__ import absolute_import, division, print_function
import scitbx.math
from scitbx import differential_evolution as de
from scitbx import simplex
from scitbx import lbfgs
from scitbx import direct_search_simulated_annealing as dssa
from scitbx.array_family import flex
from scitbx.stdlib import math
import time, sys
from fractions import Fraction
from six.moves import range
from six.moves import zip

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
  lfg =  scitbx.math.log_factorial_generator(nmax)
  #rzfa = scitbx.math.zernike_2d_radial(n,l,lfg)
  #rap = scitbx.math.zernike_2d_polynome(n,l,rzfa)
  rap = scitbx.math.zernike_2d_polynome(n,l)#,rzfa)

  image = flex.vec3_double()

  original=open('original.dat','w')
  count = 0
  for x in range(-N, N+1):
    for y in range(-N, N+1):
      rr = math.sqrt(x*x+y*y)/N
      if rr>1.0:
        value=0.0
      else:
        tt = math.atan2(y,x)
        value = rap.f(rr,tt)
        value = value.real
        count = count + 1
      image.append([x+N,y+N,value])
      print(x+N,y+N, value, file=original)
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

  NP=int(math.sqrt( image.size() ))
  N=NP/2
  grid_2d = scitbx.math.two_d_grid(N, nmax)
  grid_2d.clean_space( image )
  grid_2d.construct_space_sum()
  tt2=time.time()
  print("time used: ", tt2-tt1)
  zernike_2d_mom = scitbx.math.two_d_zernike_moments( grid_2d, nmax )

  moments = zernike_2d_mom.moments()

  tt2=time.time()
  print("time used: ", tt2-tt1)

  coefs = flex.real( moments )
  nl_array = scitbx.math.nl_array( nmax )
  nls = nl_array.nl()
  nl_array.load_coefs( nls, coefs )
  lfg =  scitbx.math.log_factorial_generator(nmax)

  print(nl_array.get_coef(n,l)*2)

  for nl, c in zip( nls, moments):
    if(abs(c)<1e-3):
      c=0
    print(nl, c)

  print()
  reconst=flex.complex_double(NP**2, 0)
  for nl,c in zip( nls, moments):
    n=nl[0]
    l=nl[1]
    if(l>0):
      c=c*2
    #rzfa = scitbx.math.zernike_2d_radial(n,l,lfg)
    rap = scitbx.math.zernike_2d_polynome(n,l) #,rzfa)
    i=0
    for x in range(0,NP):
      x=x-N
      for y in range(0,NP):
        y=y-N
        rr = math.sqrt(x*x+y*y)/N
        if rr>1.0:
          value=0.0
        else:
          tt = math.atan2(y,x)
          value = rap.f(rr,tt)
        reconst[i]=reconst[i]+value*c
        i=i+1

  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      if(value>0):
        print(x,y,image[i][2],value, file=rebuilt)
      i=i+1


  rebuilt.close()


def tst_2d_poly(n,l):
  nmax=max(n,20)
  np=50
  x,y=0.1,0.9
  r,t=math.sqrt(x*x+y*y),math.atan2(y,x)
  lfg =  scitbx.math.log_factorial_generator(nmax)
  #rzfa = scitbx.math.zernike_2d_radial(n,l,lfg)
  #rap = scitbx.math.zernike_2d_polynome(n,l,rzfa)
  rap = scitbx.math.zernike_2d_polynome(n,l)#,rzfa)
  rt_value=rap.f(r,t)
  grid = scitbx.math.two_d_grid(np, nmax)
  zm2d = scitbx.math.two_d_zernike_moments(grid, nmax)
  xy_value=zm2d.zernike_poly(n,l,x,y)

  print(rt_value, xy_value, abs(rt_value), abs(xy_value))

def tst_2d_zm(n,l):
  nmax=max(n,20)
  np=100
  points=flex.double(range(-np,np+1))/np
  grid = scitbx.math.two_d_grid(np, nmax)
  zm2d = scitbx.math.two_d_zernike_moments(grid, nmax)
  image = flex.vec3_double()

  output=file('testmap.dat','w')

  for x in points:
    for y in points:
      r=math.sqrt(x*x+y*y)
      if(r>1.0):
        value=0.0
      else:
        value=zm2d.zernike_poly(n,l,x,y).real
      image.append([x*np+np,y*np+np, value])

  grid.clean_space( image )
  grid.construct_space_sum()
  zernike_2d_mom = scitbx.math.two_d_zernike_moments( grid, nmax )

  moments = zernike_2d_mom.moments()

  coefs = flex.real( moments )
  nl_array = scitbx.math.nl_array( nmax )
  nls = nl_array.nl()

  for nl, c in zip( nls, moments):
    if(abs(c)<1e-3):
      c=0
    print(nl, c)

  NP=np*2+1
  reconst=zernike_2d_mom.zernike_map(nmax, np)
  i = 0
  for x in range(0,NP):
    for y in range(0,NP):
      value=reconst[i].real
      if(value>0):
        print(x,y,image[i][2],value, file=output)
      i=i+1


  output.close()


def integrate_triple_zernike2d(n1,n2,n3,m, Bnmk_obj):
  value=Fraction(0)
  temp = long(0)
  ck = [long(0)]*(n1+n2+n3+1)
  for k1 in range(m,n1+1,2):
    for k2 in range(m,n2+1,2):
      for k3 in range(m,n3+1,2):
       # value = value+Bnmk_obj.get_coef(n1,m,k1)*Bnmk_obj.get_coef(n2,m,k2)*Bnmk_obj.get_coef(n3,m,k3)/(k1+k2+k3+2.0)
        temp = Bnmk_obj.get_coef(n1,m,k1)*Bnmk_obj.get_coef(n2,m,k2)*Bnmk_obj.get_coef(n3,m,k3)
        ck[k1+k2+k3] = ck[k1+k2+k3]+temp

  for kk in range(3*m,n1+n2+n3+1,2):
  #  print "%4d, %30d"%(kk,ck[kk])
    value = value + Fraction( ck[kk],(kk+2))
  return float(value)

class Bnmk(object):
  "Bnmk coefficient object hold 2d zernike expansion coefs"
  def __init__(self, nmax):
    self.nmax=nmax
    self.Bnmk=scitbx.math.nmk_array(nmax)
    self.initialize_bnmk()

  def initialize_bnmk(self):
    for n in range(self.nmax, -1,-1):
      self.Bnmk.set_coef(n,n,n,1.0)
      for m in range(n-2,-1,-2):
        value = self.Bnmk.get_coef(n,m+2,n)*(n+m+2.0)/(n-m)
        self.Bnmk.set_coef(n,m,n,value)
        for k in range(n-2,m-1,-2):
          value = -self.Bnmk.get_coef(n,m,k+2)*(k+m+2.0)*(k+2.0-m)/(k+2.0+n)/(n-k)
          self.Bnmk.set_coef(n,m,k,value)

  def get_coef(self,n,m,k):
    return int(self.Bnmk.get_coef(n,m,k).real)

  def print_bnmk(self):
    for n in range(self.nmax+1):
      for m in range(n,-1,-2):
        for k in range(m,n+1,2):
          print(n,m,k,self.get_coef(n,m,k))


def tst_integral_triple_zernike2d(nmax):
  Bnmk_obj = Bnmk(nmax)
  #Bnmk_obj.print_bnmk()
  coef_table = []

  for m in range(nmax+1):
    C_m_3n = scitbx.math.nmk_array(nmax)
    for n1 in range(m,nmax+1,2):
      for n2 in range(m,n1+1,2):
        for n3 in range(m,n2+1,2):
          value = integrate_triple_zernike2d(n1,n2,n3,m,Bnmk_obj)
          C_m_3n.set_coef(n1,n2,n3,value)
        #  print m,n1,n2,n3,value.real
    coef_table.append( C_m_3n )
  return coef_table

def calc_Cnm_4m_from_Inm(m, nmax, Inm, coef_table_m):
  Cnm_m=flex.double()
  for n in range(m, nmax+1, 2 ):
    temp = 0
    for n1 in range( m,nmax+1,2 ):
      n1_indx = (n1-m)/2
      for n2 in range( m,nmax+1,2 ):
        n2_indx = (n2-m)/2
        i,j,k=sorted([n,n1,n2],reverse=True)
        temp = temp + coef_table_m.get_coef(i,j,k).real*Inm[n1_indx]*Inm[n2_indx]
    Cnm_m.append( temp )
  return Cnm_m


def calc_Cnm_from_Inm( Inm, coef_table, nmax ):
  Cnm=scitbx.math.nl_array(nmax)
  for n in range( 0,nmax+1,2 ): # only even number (n,m)
    for m in range(0,n+1,2 ):
      temp = 0
      for n1 in range( m,nmax+1,2 ):
        for n2 in range( m,n1,2 ):
          i,j,k=sorted([n,n1,n2],reverse=True)
          temp = temp + coef_table[m].get_coef(i,j,k).real*Inm.get_coef(n1,m)*Inm.get_coef(n2,m)

        i,j,k=sorted([n,n1,n1],reverse=True)
        temp = temp + coef_table[m].get_coef(i,j,k).real*Inm.get_coef(n1,m)**2.0/2.0
      Cnm.set_coef(n,m,temp*2.0)
  return Cnm


def comp_Cnm_calculations(nmax):
  Inm=scitbx.math.nl_array(nmax)
  nls = Inm.nl()
  size = nls.size()
  coefs = flex.random_double(size)
  Inm.load_coefs( nls, coefs )
  coef_table = tst_integral_triple_zernike2d(nmax)
  Cnm = calc_Cnm_from_Inm(Inm, coef_table, nmax )

  for mm in range(0,nmax+1,2):
    coef_table_m = coef_table[mm]
    Inm_m=flex.double()
    Cnm_m=flex.double()
    for nn in range(mm,nmax+1,2):
      Inm_m.append( Inm.get_coef(nn,mm) )
      Cnm_m.append( Cnm.get_coef(nn,mm) )
    Cnm_m_new = calc_Cnm_4m_from_Inm(mm,nmax,Inm_m,coef_table_m)

    for ii,jj in zip(Cnm_m, Cnm_m_new):
      assert(abs(ii-jj)<1e-6)




def tst_solving_Inm(nmax):
  Inm=scitbx.math.nl_array(nmax)
  nls = Inm.nl()
  size = nls.size()
  coefs = flex.random_double(size)
  Inm.load_coefs( nls, coefs )

  coef_table = tst_integral_triple_zernike2d(nmax)
  Cnm = calc_Cnm_from_Inm( Inm, coef_table, nmax )

  new_Inm=scitbx.math.nl_array(nmax)

  for n in range( nmax+1 ):
    for m in range(n,-1,-2):
      Cnm_value = Cnm.get_coef( n, m )
      coef = coef_table[m].get_coef(n,n,n).real
#      value = math.sqrt( Cnm_value/coef )
      if(coef != 0):
        value = ( Cnm_value/coef )
        print(n,m,Inm.get_coef(n,m)**2, value)


class inm_refine(object):
  def __init__(self, nmax, Cnm, coef_table,m):
    self.nmax=nmax
    self.Cnm = Cnm
    self.coef_table = coef_table
    self.Inm = scitbx.math.nl_array(nmax)

    self.m=m ## testing the most populated cases
    self.coef_table_m = self.coef_table[self.m]
    self.x = None
    self.ndim = (self.nmax-self.m)/2+1
    self.n = self.ndim
    self.domain =[(-1.0,1.0)]*self.ndim
    self.target_data=flex.double()
    self.n_list=range(self.m,nmax+1,2)
    self.n_indx_list=range(self.ndim)
    for nn in self.n_list:
      self.target_data.append( self.Cnm.get_coef(nn,self.m) )
#    self.solution=flex.random_double(self.n)*2-1
#    self.target_data=self.calc_Cnm_from_Inm(self.solution)
    self.scale = self.target_data.norm()**2
#    self.optimizer = de.differential_evolution_optimizer(self, population_size=self.ndim,eps=1e-8,n_cross=self.n/2)
#    self.run_simplex()
#    self.run_dssa()
    lowest_score=1e8
    for ii in range(10):
      self.run_lbfgs()
      score=self.target(self.x)
      if(score < lowest_score):
        lowest_score=score
        self.best_solution=self.x.deep_copy()
      print(list(self.x), score)
    print(lowest_score)

  def run_lbfgs(self):
    self.h = 0.000001
    self.x = flex.random_double(self.n)*2-1
    lbfgs.run(target_evaluator=self)

  def compute_functional_and_gradients(self):
   # t, dd = self.derivs( self.x , h=self.h)
    t=self.target(self.x)
    dd = self.analytic_f_prime(self.x)
    return t, dd

  def derivs(self,vector,h=0.0001):
    t  = self.target(vector)
    dd = []
    for ii in range(vector.size()):
      vector_1 = vector.deep_copy()
      vector_1[ii]=vector_1[ii]+h
      th = self.target(vector_1)
      delta = (th-t)/h
      dd.append( delta )
    delta = flex.double(dd)
    return t,delta

  def callback_after_step(self, minimizer):
    return

  def run_dssa(self):
    start_matrix=[]
    for ii in range(self.n):
      start_matrix.append( flex.random_double(self.n)*2.0-1.0 )
    start_matrix.append( flex.double(self.n,0) )
    self.optimizer = dssa.dssa( dimension=self.n, matrix=start_matrix, evaluator=self, tolerance=1e-5, further_opt=True )
    self.x = self.optimizer.get_solution()


  def run_simplex(self):
    start_matrix=[]
    for ii in range(self.n):
      start_matrix.append( flex.random_double(self.n)*2.0-1.0 )
    start_matrix.append( flex.double(self.n,0) )
    self.optimizer = simplex.simplex_opt( dimension=self.n, matrix=start_matrix, evaluator=self, tolerance=1e-5 )
    self.x = self.optimizer.get_solution()

  def analytic_f_prime(self,Inm):
    f_prime=flex.double()
    pre_factor=(self.calc_data-self.target_data)*(4.0)/self.scale

    for n in self.n_list:
      temp = 0
      for n1,n1_indx in zip(self.n_list,self.n_indx_list):
        temp1 = 0
        for n2,n2_indx in zip(self.n_list,self.n_indx_list):
          i,j,k=sorted([n,n1,n2],reverse=True)
          temp1 = temp1 + self.coef_table_m.get_coef(i,j,k).real*Inm[n2_indx]
        temp = temp+temp1*pre_factor[n1_indx]
      f_prime.append(temp)
    return f_prime


  def target(self, x):
    self.calc_data=self.calc_Cnm_from_Inm(x)
    score = flex.sum_sq( self.calc_data-self.target_data)
    return score/self.scale

  def calc_Cnm_from_Inm(self, Inm):
    Cnm=flex.double()
    for n in range(self.m,self.nmax+1,2):
      temp = 0
      for n1 in range(self.m,self.nmax+1,2):
        n1_indx = (n1-self.m)/2
        for n2 in range(self.m,self.nmax+1,2):
          n2_indx = (n2-self.m)/2
          i,j,k=sorted([n,n1,n2],reverse=True)
          temp = temp + self.coef_table_m.get_coef(i,j,k).real*Inm[n1_indx]*Inm[n2_indx]
      Cnm.append( temp )
    return Cnm

def test_solving(nmax, m):
  Inm=scitbx.math.nl_array(nmax)
  nls = Inm.nl()
  size = nls.size()
  coefs = flex.random_double(size)
  Inm.load_coefs( nls, coefs )

  coef_table = tst_integral_triple_zernike2d(nmax)
  Cnm = calc_Cnm_from_Inm( Inm, coef_table, nmax )
  t1 = time.time()
  refine_de_obj = inm_refine( nmax, Cnm, coef_table, m )
  solution=flex.double()
  for nn in range(m,nmax+1,2):
    solution.append( Inm.get_coef(nn,m) )
  for ss,xx in zip(solution, refine_de_obj.best_solution):
    print(ss,xx)
  print("score=",refine_de_obj.target(refine_de_obj.x))
  print("score=",refine_de_obj.target(solution))
  t2 = time.time()
  print("nmax=",nmax, "time used:", t2-t1)




if __name__ == "__main__":
  args = sys.argv[1:]
  if(len(args) == 4):
    m=int(args[0])
    n1=int(args[1])
    n2=int(args[2])
    n3=int(args[3])
    Bnmk_obj = Bnmk(max(n1,n2,n3))
    print(integrate_triple_zernike2d(n1,n2,n3,m,Bnmk_obj))
    exit()

  if(len(args) == 1):
    nmax=int(args[0])
  elif(len(args)==0):
    nmax=5
  if(len(args) == 2):
    nmax=int(args[0])
    m=int(args[1])
    assert m<=nmax
    test_solving(nmax,m)
  exit()
  tst_integral_triple_zernike2d(nmax)
  comp_Cnm_calculations(nmax)


  exit()


  t1=time.time()
  tst_2d_zm(41,1)
  t2=time.time()
  print("xy:  ", t2-t1)
  tst_2d_zernike_mom(41,1)
  t3=time.time()
  print("rt:  ", t3-t2)
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
  print("time used: ", t2-t1)
  print("OK")
