from __future__ import division
from libtbx import adopt_init_args
from scitbx import matrix
import math

# transformation for derivative of a unit vector
def G_unitvector(Gu, r):
  x, y, z = r[0], r[1], r[2]
  Gxu, Gyu, Gzu = Gu[0], Gu[1], Gu[2]
  #denom = math.sqrt((x**2 + y**2 + z**2)**3)
  denom = (r.length())**3
  if(denom == 0):
    raise RuntimeError("division by zero in G_unitvector")
  Gx = ( (y**2 + z**2)*Gu[0] - x*y*Gu[1] - x*z*Gu[2] ) / denom
  Gy = ( -x*y*Gu[0] + (x**2 + z**2)*Gu[1] - y*z*Gu[2] ) / denom
  Gz = ( -x*z*Gu[0] - y*z*Gu[1] + (x**2 + y**2)*Gu[2] ) / denom
  return matrix.col([Gx, Gy, Gz])
#
# transformation for cross product
def G_crossproduct(Gv, r1, r2):
  x1, y1, z1 = r1
  x2, y2, z2 = r2
  Gx_v, Gy_v, Gz_v = Gv
  Gx1 = -z2 * Gy_v + y2 * Gz_v
  Gy1 =  z2 * Gx_v - x2 * Gz_v
  Gz1 = -y2 * Gx_v + x2 * Gy_v
  Gx2 =  z1 * Gy_v - y1 * Gz_v
  Gy2 = -z1 * Gx_v + x1 * Gz_v
  Gz2 =  y1 * Gx_v - x1 * Gy_v
  Gr1 = matrix.col([Gx1, Gy1, Gz1])
  Gr2 = matrix.col([Gx2, Gy2, Gz2])
  return Gr1, Gr2
#

class modify_gradients(object):

  def __init__(self,
               sites_cart,
               h_parameterization,
               grads):
    adopt_init_args(self, locals())

    for hp in self.h_parameterization:
      if (hp == None): continue
      ih = hp.ih
      #hp = self.h_parameterization[ih]
      if(hp.htype == 'unk'):
        continue
      a0, a1, a2 = hp.a0, hp.a1, hp.a2
      a, b, h  = hp.a, hp.b, hp.h
      dh = hp.dist_h
      rh = matrix.col(sites_cart[ih])
      r0 = matrix.col(sites_cart[a0])
      GH = matrix.col(self.grads[ih])
      G0 = matrix.col(self.grads[a0])
      G1 = matrix.col(self.grads[a1])
      G2 = matrix.col(self.grads[a2])
      if (hp.htype in ['flat_2neigbs', '2neigbs', '2tetra']):
        r1 = matrix.col(sites_cart[a1])
        r2 = matrix.col(sites_cart[a2])
        r10 = r1 - r0
        r20 = r2 - r0
        u10 = (r1 - r0).normalize()
        u20 = (r2 - r0).normalize()
      # step 1
        GuH0 = dh * GH
        G01 = GH
      # step 2
        if(hp.htype == 'flat_2neigbs'):
          rh0 = (a * u10 + b * u20)
          # alternatively:
          #length = math.sqrt(rh0.dot(rh0))
          length = math.sqrt(a*a + b*b + 2*a*b*(u10).dot(u20))
          Gu10 = a/length*GuH0-a*b*u20*(1/length)**3*rh0.dot(GuH0)
          Gu20 = b/length*GuH0-a*b*u10*(1/length)**3*rh0.dot(GuH0)
        elif (hp.htype == '2tetra'):
          alpha = hp.h
          v = u10.cross(u20)
          d = a * u10 + b * u20
        # step 2
          Gd0 = math.cos(alpha) * GuH0
          Gv0 = math.sin(alpha) * GuH0
        # step 3
          Gd = G_unitvector(Gu = Gd0, r  = d)
          Gv = G_unitvector(Gu = Gv0, r  = v)
        # step 4
          Gu10_1 = a * Gd
          Gu20_1 = b * Gd
        # step 5
          Gu10_2, Gu20_2 = G_crossproduct(Gv, u10, u20)
          Gu10 = Gu10_1 + Gu10_2
          Gu20 = Gu20_1 + Gu20_2
        elif (hp.htype == '2neigbs'):
          v = u10.cross(u20)
          vl = v.length()
          v0 = v.normalize()
          rh0 = (a * u10 + b * u20 + h * v0)
          length = math.sqrt(rh0.dot(rh0))
          # alternatively:
          #length = math.sqrt(a*a + b*b + h*h + 2*a*b*(u10).dot(u20))
          Gu10_1 = a/length*GuH0 - (a*b*u20)*(1/length)**3*rh0.dot(GuH0)
          Gu20_1 = b/length*GuH0 - (a*b*u10)*(1/length)**3*rh0.dot(GuH0)
          Gv0 = h/length*GuH0
      # step 3
          Gv = G_unitvector(Gu = Gv0, r  = v)
      # step 4
          Gu10_2, Gu20_2 = G_crossproduct(Gv, u10, u20)
          Gu10 = Gu10_1 + Gu10_2
          Gu20 = Gu20_1 + Gu20_2
        # step 5
        Gr10 = G_unitvector(Gu = Gu10, r  = r10)
        Gr20 = G_unitvector(Gu = Gu20, r  = r20)
        #
        self.grads[a0] = (G0[0] + GH[0] - Gr10[0] - Gr20[0],
                          G0[1] + GH[1] - Gr10[1] - Gr20[1],
                          G0[2] + GH[2] - Gr10[2] - Gr20[2])
        self.grads[a1] = (G1[0] + Gr10[0],
                          G1[1] + Gr10[1],
                          G1[2] + Gr10[2])
        self.grads[a2] = (G2[0] + Gr20[0],
                          G2[1] + Gr20[1],
                          G2[2] + Gr20[2])
      #
      if (hp.htype == '3neigbs'):
        a3 = hp.a3
        G3 = matrix.col(self.grads[a3])
        r1 = matrix.col(sites_cart[a1])
        r2 = matrix.col(sites_cart[a2])
        r3 = matrix.col(sites_cart[a3])
        r10, r20, r30 = r1 - r0, r2 - r0, r3 - r0
        u10 = r10.normalize()
        u20 = r20.normalize()
        u30 = r30.normalize()
        rh0 = (a * u10 + b * u20 + h * u30)
        # alternatively:
        #length = math.sqrt(rh0.dot(rh0))
        length = math.sqrt(a*a + b*b + h*h+ 2*a*b*(u10).dot(u20)
          + 2*a*h*(u10).dot(u30) + 2*b*h*(u20).dot(u30))
      # step 1
        GuH0 = dh * GH
        G01 = GH
      # step 2
        Gu10 = a/length*GuH0 - (a*b*u20+a*h*u30)*(1/length)**3*rh0.dot(GuH0)
        Gu20 = b/length*GuH0 - (a*b*u10+b*h*u30)*(1/length)**3*rh0.dot(GuH0)
        Gu30 = h/length*GuH0 - (a*h*u10+b*h*u20)*(1/length)**3*rh0.dot(GuH0)
      # step 3
        Gr10 = G_unitvector(Gu = Gu10, r  = r10)
        Gr20 = G_unitvector(Gu = Gu20, r  = r20)
        Gr30 = G_unitvector(Gu = Gu30, r  = r30)
        #
        self.grads[a0] = (G0[0] + GH[0] - Gr10[0] - Gr20[0] - Gr30[0],
                          G0[1] + GH[1] - Gr10[1] - Gr20[1] - Gr30[1],
                          G0[2] + GH[2] - Gr10[2] - Gr20[2] - Gr30[2])
        self.grads[a1] = (G1[0] + Gr10[0],
                          G1[1] + Gr10[1],
                          G1[2] + Gr10[2])
        self.grads[a2] = (G2[0] + Gr20[0],
                          G2[1] + Gr20[1],
                          G2[2] + Gr20[2])
        self.grads[a3] = (G3[0] + Gr30[0],
                          G3[1] + Gr30[1],
                          G3[2] + Gr30[2])

      if (hp.htype in ['alg1b', 'prop', 'alg1a']):
        r1 = matrix.col(sites_cart[a1])
        rb1 = matrix.col(sites_cart[a2])
        alpha, phi, n = hp.a, hp.b, hp.n
        k1 = -math.cos(alpha)
        k2 = math.sin(alpha)*math.cos(phi + n*2*math.pi/3)
        k3 = math.sin(alpha)*math.sin(phi + n*2*math.pi/3)
        r01 = r0 - r1
        rb10 = rb1 - r1
        v1 = r0 - r1
        u1 = (r0 - r1).normalize()
        v2 = rb10 - ((rb10).dot(u1)) * u1
        u2 = v2.normalize()
        u3 = u1.cross(u2)
      # step 1
        GuH0 = dh * GH
        G01 = GH
      # step 2
        Gu1_1 = k1 * GuH0
        Gu2_1 = k2 * GuH0
        Gu3   = k3 * GuH0
      # step 3
        Gu1_2, Gu2_2 = G_crossproduct(Gu3, u1, u2)
        Gu2 = Gu2_1 + Gu2_2
      # step 4
        Gv2 = G_unitvector(Gu = Gu2, r  = v2)
      # step 5
        Gu1 = Gu1_1 + Gu1_2 - (rb10.dot(u1))*Gv2 - (u1.dot(Gv2))*rb10
        Grb1 = Gv2 - (u1.dot(Gv2))*u1
      # step 6
        Gr01 = G_unitvector(Gu = Gu1, r  = r01)
      #
        self.grads[a0] = (G0[0] + GH[0] + Gr01[0],
                          G0[1] + GH[1] + Gr01[1],
                          G0[2] + GH[2] + Gr01[2])
        self.grads[a1] = (G1[0] - Gr01[0] - Grb1[0],
                          G1[1] - Gr01[1] - Grb1[1],
                          G1[2] - Gr01[2] - Grb1[2])
        self.grads[a2] = (G2[0] + Grb1[0],
                          G2[1] + Grb1[1],
                          G2[2] + Grb1[2])
      #
      #self.grads[ih] = (0,0,0)
