from __future__ import absolute_import, division, print_function
from cctbx.array_family import flex
import math, sys, random
from libtbx.test_utils import approx_equal
from libtbx.utils import format_cpu_times
import random, copy
from cctbx import sgtbx
from cctbx import adptbx
from cctbx.development import random_structure
import scitbx.math

from mmtbx.tls import tools
from mmtbx_tls_ext import *
from six.moves import range

if (1): # fixed random seed to avoid rare failures
  random.seed(0)
  flex.set_random_seed(0)

def branch_3_mn():
  m = None
  n = None
  small = 1.e-15
  for i in range(10000):
      m2=[]
      n2=[]
      for i in range(4):
          r1 = random.random()
          r2 = random.random()
          r3 = random.random()
          if(r3 > 0.1):
             r1 = r2
          m2.append(r1)
          n2.append(r2)
      p1 = 0.5 * (m2[0]+m2[3] + math.sqrt(4*m2[1]*m2[2]+(m2[0]-m2[3])**2))
      p2 = 0.5 * (m2[0]+m2[3] - math.sqrt(4*m2[1]*m2[2]+(m2[0]-m2[3])**2))
      q1 = 0.5 * (n2[0]+m2[3] + math.sqrt(4*n2[1]*n2[2]+(n2[0]-n2[3])**2))
      q2 = 0.5 * (n2[0]+m2[3] - math.sqrt(4*n2[1]*n2[2]+(n2[0]-n2[3])**2))
      if(min(p1,p2) > 0.0 and min(q1,q2) > 0.0):
          r = random.random()
          if(r > 0.5):
             r1 = r3
             r2 = r2
          m = [m2[0],m2[3],r1,m2[1],0,0]
          n = [n2[0],n2[3],r2,n2[1],0,0]
          if([adptbx.is_positive_definite(m),
              adptbx.is_positive_definite(n)].count(True)==2):
              esm = adptbx.eigensystem(m)
              esn = adptbx.eigensystem(n)
              vn = esn.values()
              vm = esm.values()
              mmin = flex.min(flex.double(vm))
              nmin = flex.min(flex.double(vn))
              if(abs(abs(mmin) - abs(nmin)) < small and mmin> 0. and nmin> 0.):
                  for i, v in enumerate(vn):
                      if(abs(abs(nmin) - v) < small): break
                  for j, v in enumerate(vm):
                      if(abs(abs(mmin) - v) < small): break
                  vecn = flex.double(esn.vectors(i))
                  vecm = flex.double(esm.vectors(j))
                  if(flex.mean(vecm-vecn) < small): break
              else:
                m = None
                n = None
  assert [m,n] != [None,None]
  assert [adptbx.is_positive_definite(m),
          adptbx.is_positive_definite(n)].count(True)==2
  r = random.random()
  if(r > 0.5):
     m = adptbx.random_rotate_ellipsoid(u_cart = m)
     n = adptbx.random_rotate_ellipsoid(u_cart = n)
  return m,n

def branch_2_mn(small):
  m=[]
  n=[]
  while (len(n)==0 or len(m)==0):
    r1 = random.random()
    r2 = random.random()
    r3 = random.random()
    if(r1 >= r2 and r2 >= r3 and r3>  0.0):
       r = random.random()
       if(r > 0.5):
          m = [r1,r2,r3,0,0,0]
       else:
          m = [r1,r1,r1,0,0,0]
       assert adptbx.is_positive_definite(m,small)
    if(r1 >= r2 and r2 >= r3 and r3 >  0.0):
       r = random.random()
       if(r > 0.5):
          n = [r1,r3,r2,0,0,0]
       else:
          n = [r1,r1,r1,0,0,0]
       assert adptbx.is_positive_definite(n,small)
  r = random.random()
  if(r > 0.5):
     m = adptbx.random_rotate_ellipsoid(u_cart = m)
     n = adptbx.random_rotate_ellipsoid(u_cart = n)
  return m,n

def zero_mat(x):
  for i_seq, v in enumerate(x):
    if(abs(v) < 1.e-4): x[i_seq]=0.0
  return x

def factor(x, i):
  if(i < 50000): scale = 0.01
  else: scale = 10.0
  for i_seq, v in enumerate(x):
    x[i_seq]=v*scale
  return x

def exercise_branch_2_1(small = 1.e-9):
  m = [2,1,2,0,0,0]
  n = [2,1,1,0,0,0]
  counter = 0
  trials = 100000
  i = 0
  branch_0       = 0
  branch_1       = 0
  branch_1_1     = 0
  branch_1_2     = 0
  branch_1_2_1   = 0
  branch_1_2_2   = 0
  branch_1_2_3   = 0
  branch_1_2_3_1 = 0
  branch_1_2_3_2 = 0
  while i < trials:
      i += 1
      m = [2,1,2,0,0,0]
      n = [2,1,1,0,0,0]
      counter += 1
      c = scitbx.math.euler_angles_as_matrix(
                    [random.uniform(0,360) for j in range(3)], deg=True).elems
      r = random.random()
      if(r<0.25):
         m = adptbx.c_u_c_transpose(c, m)
         n = adptbx.c_u_c_transpose(c, n)
      elif(r>=0.25 and r < 0.5):
         m = adptbx.c_u_c_transpose(c, m)
      elif(r>=0.5 and r < 0.75):
         n = adptbx.c_u_c_transpose(c, n)
      else:
         r = random.random()
         run_away = 1000
         while run_away > 0:
            run_away -= 1
            r = random.random()
            if(r<0.33):
               m = adptbx.c_u_c_transpose(c, m)
               n = adptbx.c_u_c_transpose(c, n)
            elif(r>=0.33 and r < 0.66):
               m = adptbx.c_u_c_transpose(c, m)
            else:
               n = adptbx.c_u_c_transpose(c, n)
            m = adptbx.c_u_c_transpose(c, m)
            n = adptbx.c_u_c_transpose(c, n)
            m = [m[0],m[1],m[2],m[3]*random.choice((0,1)),m[4]*random.choice((0,1)),m[5]*random.choice((0,1))]
            n = [n[0],n[1],n[2],n[3]*random.choice((0,1)),n[4]*random.choice((0,1)),n[5]*random.choice((0,1))]
            m = factor(list(m), i)
            n = factor(list(n), i)
            if(adptbx.is_positive_definite(m,small) and
               adptbx.is_positive_definite(n,small)): break


      m_dc = copy.deepcopy(m)
      n_dc = copy.deepcopy(n)
      qq = tools.common(n,m,small)
      #assert approx_equal(m,m_dc)
      #assert approx_equal(n,n_dc)

      if(qq.branch_0()      ):   branch_0       += 1
      if(qq.branch_1()      ):   branch_1       += 1
      if(qq.branch_1_1()    ):   branch_1_1     += 1
      if(qq.branch_1_2()    ):   branch_1_2     += 1
      if(qq.branch_1_2_1()  ):   branch_1_2_1   += 1
      if(qq.branch_1_2_2()  ):   branch_1_2_2   += 1
      if(qq.branch_1_2_3()  ):   branch_1_2_3   += 1
      if(qq.branch_1_2_3_1()):   branch_1_2_3_1 += 1
      if(qq.branch_1_2_3_2()):   branch_1_2_3_2 += 1

      if (counter >= 10000):
         counter = 0
         print("."*30)
         print("i= ", i, "out of ", trials)
         print("branch_0       = ", branch_0)
         print("branch_1       = ", branch_1)
         print("branch_1_1     = ", branch_1_1)
         print("branch_1_2     = ", branch_1_2)
         print("branch_1_2_1   = ", branch_1_2_1)
         print("branch_1_2_2   = ", branch_1_2_2)
         print("branch_1_2_3   = ", branch_1_2_3)
         print("branch_1_2_3_1 = ", branch_1_2_3_1)
         print("branch_1_2_3_2 = ", branch_1_2_3_2)
         sys.stdout.flush()


def exercise_2(small = 1.e-9):
  m = [2,1,2,0,0,0]
  n = [2,2,1,0,0,0]
  m_dc = copy.deepcopy(m)
  n_dc = copy.deepcopy(n)
  qq = tools.common(n,m,small)
  assert approx_equal(m,m_dc)
  assert approx_equal(n,n_dc)
  assert qq.t() == (2.0, 1.0, 1.0, 0.0, 0.0, 0.0)


def exercise(small = 1.e-9):
  for symbol in ["P 1"]:
      space_group_info = sgtbx.space_group_info(symbol = symbol)
      random_xray_structure = random_structure.xray_structure(
                                       space_group_info  = space_group_info,
                                       elements          = ["N"]*10,
                                       volume_per_atom   = 50.0,
                                       random_u_iso      = False,
                                       u_iso             = adptbx.b_as_u(20.0))
      sg = random_xray_structure.space_group()
      uc = random_xray_structure.unit_cell()
      print(symbol, uc)
      print()
      sys.stdout.flush()
      counter = 0
      trials = 100000
      i = 0
      branch_0       = 0
      branch_1       = 0
      branch_1_1     = 0
      branch_1_2     = 0
      branch_1_2_1   = 0
      branch_1_2_2   = 0
      branch_1_2_3   = 0
      branch_1_2_3_1 = 0
      branch_1_2_3_2 = 0

      while i < trials:
          i += 1
          counter += 1
          r = random.random()
          if(r < 0.333):
             m = adptbx.random_u_cart(u_scale=20.*random.random(), u_min=0)
             n = adptbx.random_u_cart(u_scale=20.*random.random(), u_min=0)
             while 1:
               for ind in range(6):
                   r = random.random()
                   m = flex.double(m)
                   if(r > 0.5):
                      m[ind] = n[ind]
               m = list(m)
               if(adptbx.is_positive_definite(m,0) and
                  adptbx.is_positive_definite(n,0)): break
          elif(r>=0.333 and r<0.66):
             m,n = branch_3_mn()
          else:
             m,n = branch_2_mn(0)
          m_dc = copy.deepcopy(m)
          n_dc = copy.deepcopy(n)
          m = factor(list(m), i)
          n = factor(list(n), i)
          qq = tools.common(n,m,small)
          #assert approx_equal(m,m_dc)
          #assert approx_equal(n,n_dc)

          if(qq.branch_0()      ):   branch_0       += 1
          if(qq.branch_1()      ):   branch_1       += 1
          if(qq.branch_1_1()    ):   branch_1_1     += 1
          if(qq.branch_1_2()    ):   branch_1_2     += 1
          if(qq.branch_1_2_1()  ):   branch_1_2_1   += 1
          if(qq.branch_1_2_2()  ):   branch_1_2_2   += 1
          if(qq.branch_1_2_3()  ):   branch_1_2_3   += 1
          if(qq.branch_1_2_3_1()):   branch_1_2_3_1 += 1
          if(qq.branch_1_2_3_2()):   branch_1_2_3_2 += 1


          if (counter >= 10000):
             counter = 0
             print("."*30, symbol)
             print("i= ", i, "out of ", trials)
             print("branch_0       = ", branch_0)
             print("branch_1       = ", branch_1)
             print("branch_1_1     = ", branch_1_1)
             print("branch_1_2     = ", branch_1_2)
             print("branch_1_2_1   = ", branch_1_2_1)
             print("branch_1_2_2   = ", branch_1_2_2)
             print("branch_1_2_3   = ", branch_1_2_3)
             print("branch_1_2_3_1 = ", branch_1_2_3_1)
             print("branch_1_2_3_2 = ", branch_1_2_3_2)
             sys.stdout.flush()

def exercise_x1():
  u_iso = flex.double([1.0,2.0,3.0,4.0,5.0])
  t = tools.t_from_u_cart(u_iso, 1.e-9)
  assert t == (1.0, 1.0, 1.0, 0.0, 0.0, 0.0)
  u_iso = flex.double([7.0,2.0,3.0,9.0,5.0])
  t = tools.t_from_u_cart(u_iso, 1.e-9)
  assert t == (2.0, 2.0, 2.0, 0.0, 0.0, 0.0)

def run():
  exercise_x1()
  exercise()
  exercise_branch_2_1()
  exercise_2()
  print("OK: ",format_cpu_times())

if (__name__ == "__main__"):
  run()
