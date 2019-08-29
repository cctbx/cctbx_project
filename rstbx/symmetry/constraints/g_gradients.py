from __future__ import absolute_import, division, print_function
from six.moves import range
from math import sqrt,cos,sin
from cctbx.array_family import flex

class g_gradients:
  """Use the chain rule to convert from Aij gradients to gi gradients"""

  def __init__(self, agadaptor, symred):
    """receive instances of helper classes"""
    self.agadaptor = agadaptor
    try:
      self.symred_constraints = symred.constraints
    except Exception:
      self.symred_constraints = symred


  def get_all_da(self):
    # Aij : 9 elements of the reciprocal A* matrix [Rsymop paper, eqn(3)]
    # df_Aij: 9 gradients of the functional with respect to Aij elements

    # G=(g0,g1,g2,g3,g4,g5) 6 elements of the symmetrized metrical matrix,
    #   with current increments already applied ==(a*.a*,b*.b*,c*.c*,a*.b*,a*.c*,b*.c*)
    g = self.agadaptor.G

    # The A* matrix calculated from G (lower triangular only; not general):
    #      a*_x     0    0        astrx   0     0
    # A* = a*_y   b*_y   0     == astry bstry   0
    #      a*_z   b*_z   c*_z     astrz bstrz cstrz
    #
    # Conversions from G to A*:  cstrz = sqrt(g2)
    #                            bstrz = g5/cstrz
    #                            astrz = g4/cstrz
    #                            bstry = sqrt(g1-bstrz^2)
    #                            astry = (g3-astrz*bstrz)/bstry
    #                            astrx = sqrt(g0 - astry^2 -astrz^2)
    #
    # Conversions from A* to G:  g0 = astrx^2 + astry^2 + astrz^2
    #                            g1 = bstry^2 + bstrz^2
    #                            g2 = cstrz^2
    #                            g3 = astry*bstry + astrz*bstrz
    #                            g4 = astrz*cstrz
    #                            g5 = bstrz*cstrz

    # the angles f = phi, p = psi, t = theta, along with convenient trig
    # expressions
    f = self.agadaptor.phi
    p = self.agadaptor.psi
    t = self.agadaptor.theta

    cosf = cos(f); sinf = sin(f)
    cosp = cos(p); sinp = sin(p)
    cost = cos(t); sint = sin(t)

    trig1 =  cosf*cost*sinp-sinf*sint
    trig2 =  cost*sinf+cosf*sinp*sint
    trig3 = -cost*sinf*sinp-cosf*sint
    trig4 =  cosf*cost-sinf*sinp*sint

    cstrz = sqrt(g[2])
    bstrz = g[5]/cstrz
    astrz = g[4]/cstrz
    bstry = sqrt(g[1]-bstrz*bstrz)
    astry = (g[3]-astrz*bstrz)/bstry
    astrx = sqrt(g[0] - astry*astry -astrz*astrz)

    # dAij_dphi: Partial derivatives of the 9 A components with respect to phi
    dAij_dphi = flex.double((
      -astry * trig1 - astrx * trig2 -astrz * cosf * cosp,
      -bstry * trig1 - bstrz * cosf * cosp,
      -cstrz * cosf * cosp,
      0,
      0,
      0,
      astry * trig3 + astrx * trig4 -astrz * cosp * sinf,
      bstry * trig3 - bstrz * cosp * sinf,
      -cstrz * cosp * sinf
    ))

    # dAij_dpsi: Partial derivatives of the 9 A components with respect to psi
    dAij_dpsi = flex.double((
      -astry * cosp * cost * sinf + astrz * sinf * sinp - astrx * cosp * sinf * sint,
      -bstry*cosp*cost*sinf + bstrz*sinf*sinp,
      cstrz * sinf * sinp,
      -astrz * cosp - astry * cost * sinp - astrx * sinp * sint,
      -bstrz * cosp - bstry * cost * sinp,
      -cstrz * cosp,
      astry * cosf * cosp * cost - astrz * cosf * sinp + astrx * cosf * cosp * sint,
      bstry * cosf * cosp * cost - bstrz * cosf * sinp,
      -cstrz * cosf * sinp
    ))

    # dAij_dtheta: Partial derivatives of the 9 A components with respect to theta
    dAij_dtheta = flex.double((
      -astry*trig4 + astrx*trig3,
      -bstry*trig4,
      0,
      astrx*cosp*cost - astry*cosp*sint,
      -bstry*cosp*sint,
      0,
      -astry*trig2 + astrx*trig1,
      -bstry*trig2,
      0
    ))

    # Partial derivatives of the 9 A components with respect to metrical matrix elements
    g0,g1,g2,g3,g4,g5 = g
    rad3 = g0-((g2*g3*g3+g1*g4*g4-2*g3*g4*g5)/(g1*g2-g5*g5))
    sqrt_rad3 = sqrt(rad3)

    dAij_dg0 = flex.double((
      0.5*trig4/sqrt_rad3,
      0,
      0,
      0.5*cosp*sint/sqrt_rad3,
      0,
      0,
      0.5*trig2/sqrt_rad3,
      0,
      0,
    ))

    fac4 = g2*g3-g4*g5
    fac4_sq = fac4*fac4
    rad1 =  g1-g5*g5/g2
    rad1_three_half = sqrt(rad1*rad1*rad1)
    fac3 = g5*g5-g1*g2
    rad2 = -(g2*g3*g3 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5*g5)
    factor_dg1 = fac4*fac4/(fac3*fac3*sqrt(rad2))
    fac5 = g3-(g4*g5/g2)

    dAij_dg1 = flex.double((
     -0.5*fac5*trig3/rad1_three_half + 0.5*factor_dg1*trig4,
      0.5*trig3/sqrt(rad1),
      0,
     -0.5*fac5*cosp*cost/rad1_three_half + 0.5*factor_dg1*cosp*sint,
      0.5*cosp*cost/sqrt(rad1),
      0,
     -0.5*fac5*trig1/rad1_three_half + 0.5*factor_dg1*trig2,
      0.5*trig1/sqrt(rad1),
      0
    ))

    rat5_22 = g5/(g2*g2)
    fac1 = g5*(g3-g4*g5/g2)
    fac2 = (g1*g4-g3*g5)
    fac2sq = fac2*fac2

    dAij_dg2 = flex.double((
     -0.5*rat5_22*fac1*trig3/rad1_three_half + g4*rat5_22*trig3/sqrt(rad1) +
          0.5*fac2sq*trig4/(fac3*fac3*sqrt(rad2)) + 0.5*g4*cosp*sinf/pow(g2,1.5),
      0.5*rat5_22*(g5*trig3/sqrt(rad1)+sqrt(g2)*cosp*sinf),
     -0.5*cosp*sinf/sqrt(g2),

     -0.5*rat5_22*fac1*cosp*cost/rad1_three_half + g4*rat5_22*cosp*cost/sqrt(rad1) +
          0.5*g4*sinp/pow(g2,1.5) + 0.5*(fac2sq/fac3)*cosp*sint/(fac3*sqrt(rad2)),
      0.5*rat5_22*(g5*cosp*cost/sqrt(rad1)+sqrt(g2)*sinp),
     -0.5*sinp/sqrt(g2),

     -0.5*rat5_22*fac1*trig1/rad1_three_half + g4*rat5_22*trig1/sqrt(rad1) +
          0.5*fac2sq*trig2/(fac3*fac3*sqrt(rad2)) - 0.5*g4*cosf*cosp/pow(g2,1.5),
      0.5*rat5_22*(g5*trig1/sqrt(rad1)-sqrt(g2)*cosf*cosp),
      0.5*cosf*cosp/sqrt(g2)
    ))

    dAij_dg3 = flex.double((
      trig3/sqrt(rad1) + fac4*trig4/(fac3*sqrt(rad2)),
      0,
      0,
      cosp*cost/sqrt(rad1) + fac4*cosp*sint/(fac3*sqrt(rad2)),
      0,
      0,
      trig1/sqrt(rad1) + fac4*trig2/(fac3*sqrt(rad2)),
      0,
      0
    ))

    dAij_dg4 = flex.double((
    -g5*trig3/(g2*sqrt(rad1)) + fac2*trig4/(fac3*sqrt(rad2)) - cosp*sinf/sqrt(g2),
    0,
    0,
    -g5*cosp*cost/(g2*sqrt(rad1)) - sinp/sqrt(g2) + fac2*cosp*sint/(fac3*sqrt(rad2)),
    0,
    0,
    -g5*trig1/(g2*sqrt(rad1)) + fac2*trig2/(fac3*sqrt(rad2)) + cosf*cosp/sqrt(g2),
    0,
    0
    ))

    better_ratio = (fac2/fac3)*(fac4/fac3)

    dAij_dg5 = flex.double((
    fac1*trig3/(g2*rad1_three_half) -g4*trig3/(g2*sqrt(rad1)) +
      better_ratio*trig4/sqrt(rad2),
    -g5*trig3/(g2*sqrt(rad1)) - cosp*sinf/sqrt(g2),
    0,

    fac1*cosp*cost/(g2*rad1_three_half) - g4*cosp*cost/(g2*sqrt(rad1)) +
      better_ratio*cosp*sint/sqrt(rad2),
    -g5*cosp*cost/(g2*sqrt(rad1)) -sinp/sqrt(g2),
    0,

    fac1*trig1/(g2*rad1_three_half) - g4*trig1/(g2*sqrt(rad1)) +
      better_ratio*trig2/sqrt(rad2),
    -g5*trig1/(g2*sqrt(rad1))+cosf*cosp/sqrt(g2),
    0,
    ))

    all_da = [dAij_dphi, dAij_dpsi, dAij_dtheta,
              dAij_dg0, dAij_dg1, dAij_dg2, dAij_dg3, dAij_dg4, dAij_dg5]
    return all_da

  def dB_dp(self):
    #returns the partial derivatives of the B matrix with respect to each
    # independent parameter, as a list.

    # First get the partial derivatives with respect to the six metrical
    # matrix parameters:
    dB_dg = self.get_all_da()[3:9]
    values=[]
    Nindep = self.symred_constraints.n_independent_params()
    for n in range(Nindep):
      values.append(flex.double(9))
    # Now convert to partial derivatives with respect to the independent
    for x in range(9):  # loop over the 9 elements of B
      all_gradients = [t[x] for t in dB_dg]
      g_indep = self.symred_constraints.independent_gradients(
              all_gradients=tuple(all_gradients))
      for a in range(Nindep):
        values[a][x] = g_indep[a]
    return values

  def df_dgi(self,df_dAij):
    self.all_da = self.get_all_da()
    #df_dgi, Use chain rule to calculate total derivative of functional
    # with respect to the 3 angles and 6 metrical matrix elements
    df_dgi_full = flex.double((0,0,0,0,0,0,0,0,0))
    for x in range(9):
      for y in range(9):
        df_dgi_full[x] += df_dAij[y] * self.all_da[x][y]

    g_indep = self.symred_constraints.independent_gradients(
              all_gradients=tuple(df_dgi_full[3:9]))

    value = flex.double()
    for i in range(3): value.append(df_dgi_full[i])
    for i in range(len(g_indep)): value.append(g_indep[i])
    return value

"""Input file to Mathematica for calculating the above...
Fback = {{astrx,0,0},{astry,bstry,0},{astrz,bstrz,cstrz}}
B={{Cos[f],0,Sin[f]},{0,1,0},{-Sin[f],0,Cos[f]}}
BI = Simplify[Inverse[B]]
CC = {{1,0,0},{0,Cos[p],Sin[p]},{0,-Sin[p],Cos[p]}}
CCI = Simplify[Inverse[CC]]
Simplify[CC.CCI]
DD = {{Cos[t],Sin[t],0},{-Sin[t],Cos[t],0},{0,0,1}}
DDI = Simplify[Inverse[DD]]
ROTI = BI.CCI.DDI
Aback = ROTI.Fback
derivf = D[Aback,f]//. -Cos[f]Cos[t]Sin[p]+Sin[f]Sin[t] -> negtrig1 //. -Cos[t]Sin[f]-Cos[f]Sin[p]Sin[t] -> negtrig2 //. -Cos[t]Sin[f]Sin[p]-Cos[f]Sin[t] -> trig3 //. Cos[f]Cos[t]-Sin[f]Sin[p]Sin[t] -> trig4
derivp = D[Aback,p]
derivt = D[Aback,t]//. Cos[f]Cos[t]Sin[p]-Sin[f]Sin[t] -> trig1 //. -Cos[t]Sin[f]-Cos[f]Sin[p]Sin[t] -> negtrig2 //. -Cos[t]Sin[f]Sin[p]-Cos[f]Sin[t] -> trig3 //. -Cos[f]Cos[t]+Sin[f]Sin[p]Sin[t] -> negtrig4
Aback4G = ROTI.Fback //. Cos[f]Cos[t]Sin[p]-Sin[f]Sin[t] -> trig1 //. Cos[t]Sin[f]+Cos[f]Sin[p]Sin[t] -> trig2 //. -Cos[t]Sin[f]Sin[p]-Cos[f]Sin[t] -> trig3 //. Cos[f]Cos[t]-Sin[f]Sin[p]Sin[t] -> trig4
cstrz = Sqrt[g2]
bstrz = g5/cstrz
astrz = g4/cstrz
bstry = Sqrt[g1-bstrz^2]
astry = (g3 - astrz*bstrz)/bstry
astrx = Sqrt[g0 - astry^2-astrz^2]
Simplify[D[Aback4G,g2]]//. g1-g5^2/g2 -> rad1 //. (g5^2-g1*g2)->fac3 //. g5^2*(g3-g4*g5/g2) ->g5fac1 //. -(g2*g3^2 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5^2) -> rad2 //. (g1*g4-g3*g5)^2->fac2sq
Simplify[D[Aback4G,g5]]//. g1-g5^2/g2 -> rad1 //. (g5^2-g1*g2)->fac3  //. g5*(g3-g4*g5/g2) ->fac1 //. -(g2*g3^2 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5^2) -> rad2 //. (g1*g4-g3*g5)->fac2 //. g2*g3-g4*g5 -> fac4
Simplify[D[Aback4G,g0]]//. g1-g5^2/g2 -> rad1 //. (g5^2-g1*g2)->fac3  //. g5*(g3-g4*g5/g2) ->fac1 //. -(g2*g3^2 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5^2) -> rad2 //. (g1*g4-g3*g5)->fac2 //. g2*g3-g4*g5 -> fac4 //. g0-((g2*g3^2+g1*g4^2-2*g3*g4*g5)/(g1*g2-g5^2))->rad3
Simplify[D[Aback4G,g1]]//. g1-g5^2/g2 -> rad1 //. (g5^2-g1*g2)->fac3  //. g5*(g3-g4*g5/g2) ->fac1 //. -(g2*g3^2 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5^2) -> rad2 //. (g1*g4-g3*g5)->fac2 //. g2*g3-g4*g5 -> fac4 //. g3-(g4*g5/g2) -> fac5
Simplify[D[Aback4G,g3]]//. g1-g5^2/g2 -> rad1 //. (g5^2-g1*g2)->fac3  //. g5*(g3-g4*g5/g2) ->fac1 //. -(g2*g3^2 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5^2) -> rad2 //. (g1*g4-g3*g5)->fac2 //. g2*g3-g4*g5 -> fac4
Simplify[D[Aback4G,g4]]//. g1-g5^2/g2 -> rad1 //. (g5^2-g1*g2)->fac3  //. g5*(g3-g4*g5/g2) ->fac1 //. -(g2*g3^2 +g4*(g1*g4-2*g3*g5)+g0*fac3)/(g1*g2-g5^2) -> rad2 //. (g1*g4-g3*g5)->fac2 //. g2*g3-g4*g5 -> fac4
"""
def finite_difference_test(orient):
  from rstbx.symmetry.constraints import AGconvert as AG
  from labelit.symmetry.metricsym.a_g_conversion import pp
  from libtbx.tst_utils import approx_equal
  adaptor = AG()
  adaptor.forward(orient)
  grad = g_gradients(adaptor,symred=None)
  epsilon = 1.E-10

  dAij_dphi = grad.get_all_da()[0]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi+x*epsilon,adaptor.psi,adaptor.theta)
    AGback.G = adaptor.G
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dphi: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dphi))+"\n"+\
         pp(diff_mat)
  if not ( approx_equal(dAij_dphi,diff_mat,1.E-7)): raise Exception(rule)

  dAij_dpsi = grad.get_all_da()[1]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi+x*epsilon,adaptor.theta)
    AGback.G = adaptor.G
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dpsi: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dpsi))+"\n"+\
         pp(diff_mat)
  if not( approx_equal(dAij_dpsi,diff_mat,1.E-7)): raise Exception(rule)

  dAij_dtheta = grad.get_all_da()[2]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta+x*epsilon)
    AGback.G = adaptor.G
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dtheta: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dtheta))+"\n"+\
         pp(diff_mat)
  if not( approx_equal(dAij_dtheta,diff_mat,1.E-7)): raise Exception(rule)

  g0,g1,g2,g3,g4,g5 = adaptor.G

  dAij_dg0 = grad.get_all_da()[3]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta)
    AGback.G = (g0+x*epsilon,g1,g2,g3,g4,g5)
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dg0: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dg0))+"\n"+\
         pp(diff_mat)
  if not ( approx_equal(dAij_dg0,diff_mat,1.E-7)): raise Exception(rule)

  dAij_dg1 = grad.get_all_da()[4]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta)
    AGback.G = (g0,g1+x*epsilon,g2,g3,g4,g5)
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dg1: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dg1))+"\n"+\
         pp(diff_mat)
  if not ( approx_equal(dAij_dg1,diff_mat,1.E-7)): raise Exception(rule)

  dAij_dg2 = grad.get_all_da()[5]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta)
    AGback.G = (g0,g1,g2+x*epsilon,g3,g4,g5)
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dg2: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dg2))+"\n"+\
         pp(diff_mat)
  if not( approx_equal(dAij_dg2,diff_mat,1.E-6)): raise Exception(rule)

  dAij_dg3 = grad.get_all_da()[6]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta)
    AGback.G = (g0,g1,g2,g3+x*epsilon,g4,g5)
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dg3: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dg3))+"\n"+\
         pp(diff_mat)
  if not( approx_equal(dAij_dg3,diff_mat,1.E-7)): raise Exception(rule)

  dAij_dg4 = grad.get_all_da()[7]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta)
    AGback.G = (g0,g1,g2,g3,g4+x*epsilon,g5)
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dg4: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dg4))+"\n"+\
         pp(diff_mat)
  if not( approx_equal(dAij_dg4,diff_mat,1.E-7)): raise Exception(rule)

  dAij_dg5 = grad.get_all_da()[8]
  AGback = AG()
  F = []
  for x in [-1.,1.]:
    AGback.setAngles(adaptor.phi,adaptor.psi,adaptor.theta)
    AGback.G = (g0,g1,g2,g3,g4,g5+x*epsilon)
    F.append(flex.double(AGback.back()))
  diff_mat = (F[1]-F[0])/(2.*epsilon)
  rule = "dAij_dg5: Analytical gradient vs. finite difference gradient\n"+\
         pp(list(dAij_dg5))+"\n"+\
         pp(diff_mat)
  if not( approx_equal(dAij_dg5,diff_mat,1.E-6)): raise Exception(rule)
