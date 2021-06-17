from __future__ import absolute_import, division, print_function
import scitbx.array_family.flex # import dependency

import boost_adaptbx.boost.python as bp
ext = bp.import_ext("scitbx_rigid_body_ext")
from scitbx_rigid_body_ext import *
import math
from scitbx import matrix

class pickle_import_trigger(object): pass

@bp.inject_into(ext.tardy_model)
class _():

  def __getinitargs__(O):
    return (
      O.labels,
      O.sites,
      O.masses,
      O.tardy_tree,
      O.potential_obj,
      O.near_singular_hinges_angular_tolerance_deg)

  def __getstate__(O):
    version = 2
    return (version, pickle_import_trigger(), O.pack_q(), O.pack_qd())

  def __setstate__(O, state):
    assert len(state) >= 3
    version = state[0]
    if   (version == 1): assert len(state) == 3
    elif (version == 2): assert len(state) == 4
    else: raise RuntimeError("Unknown version of pickled state.")
    O.unpack_q(q_packed=state[-2])
    O.unpack_qd(qd_packed=state[-1])

  def minimization(O, max_iterations=None, callback_after_step=None):
    from scitbx.rigid_body.essence.tardy import refinery
    return refinery(
      tardy_model=O,
      max_iterations=max_iterations,
      callback_after_step=callback_after_step)

def euler(phi, psi, the, convention):
  if(convention == "zyz"):
     result = rb_mat_zyz(the=the, psi=psi, phi=phi)
  elif(convention == "xyz"):
     result = rb_mat_xyz(the=the, psi=psi, phi=phi)
  else:
     raise RuntimeError("\nWrong rotation convention\n")
  return result

class rb_mat_zyz(object):

   def __init__(self, the, psi, phi):
     the = the * math.pi/180
     psi = psi * math.pi/180
     phi = phi * math.pi/180
     self.c_psi = math.cos(psi)
     self.c_the = math.cos(the)
     self.c_phi = math.cos(phi)
     self.s_psi = math.sin(psi)
     self.s_the = math.sin(the)
     self.s_phi = math.sin(phi)

   def rot_mat(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 =  c_the*c_psi*c_phi - s_the*s_phi
     r12 = -c_the*c_psi*s_phi - s_the*c_phi
     r13 =  c_the*s_psi
     r21 =  s_the*c_psi*c_phi + c_the*s_phi
     r22 = -s_the*c_psi*s_phi + c_the*c_phi
     r23 =  s_the*s_psi
     r31 = -s_psi*c_phi
     r32 =  s_psi*s_phi
     r33 =  c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_the(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -s_the*c_psi*c_phi - c_the*s_phi
     r12 =  s_the*c_psi*s_phi - c_the*c_phi
     r13 = -s_the*s_psi
     r21 =  c_the*c_psi*c_phi - s_the*s_phi
     r22 = -c_the*c_psi*s_phi - s_the*c_phi
     r23 =  c_the*s_psi
     r31 = 0.0
     r32 = 0.0
     r33 = 0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_psi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -c_the*s_psi*c_phi
     r12 =  c_the*s_psi*s_phi
     r13 =  c_the*c_psi
     r21 = -s_the*s_psi*c_phi
     r22 =  s_the*s_psi*s_phi
     r23 =  s_the*c_psi
     r31 = -c_psi*c_phi
     r32 =  c_psi*s_phi
     r33 = -s_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_phi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -c_the*c_psi*s_phi - s_the*c_phi
     r12 = -c_the*c_psi*c_phi + s_the*s_phi
     r13 =  0.0
     r21 = -s_the*c_psi*s_phi + c_the*c_phi
     r22 = -s_the*c_psi*c_phi - c_the*s_phi
     r23 =  0.0
     r31 =  s_psi*s_phi
     r32 =  s_psi*c_phi
     r33 =  0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

class rb_mat_xyz(object):
   """ R = Rx(the)Ry(psi)Rz(phi) """

   def __init__(self, phi, psi, the, deg=True):
     if deg:
       phi = phi * math.pi/180
       psi = psi * math.pi/180
       the = the * math.pi/180
     self.c_psi = math.cos(psi)
     self.c_phi = math.cos(phi)
     self.c_the = math.cos(the)
     self.s_psi = math.sin(psi)
     self.s_phi = math.sin(phi)
     self.s_the = math.sin(the)

   def rot_mat(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 =  c_psi*c_phi
     r12 = -c_psi*s_phi
     r13 =  s_psi
     r21 =  c_the*s_phi + s_the*s_psi*c_phi
     r22 =  c_the*c_phi - s_the*s_psi*s_phi
     r23 = -s_the*c_psi
     r31 =  s_the*s_phi - c_the*s_psi*c_phi
     r32 =  s_the*c_phi + c_the*s_psi*s_phi
     r33 =  c_the*c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_the(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 =  0.0
     r12 =  0.0
     r13 =  0.0
     r21 = -s_the*s_phi+c_the*s_psi*c_phi
     r22 = -s_the*c_phi-c_the*s_psi*s_phi
     r23 = -c_the*c_psi
     r31 =  c_the*s_phi+s_the*s_psi*c_phi
     r32 =  c_the*c_phi-s_the*s_psi*s_phi
     r33 = -s_the*c_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_psi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -s_psi*c_phi
     r12 =  s_psi*s_phi
     r13 =  c_psi
     r21 =  s_the*c_psi*c_phi
     r22 = -s_the*c_psi*s_phi
     r23 =  s_the*s_psi
     r31 = -c_the*c_psi*c_phi
     r32 =  c_the*c_psi*s_phi
     r33 = -c_the*s_psi
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm

   def r_phi(self):
     c_psi = self.c_psi
     c_the = self.c_the
     c_phi = self.c_phi
     s_psi = self.s_psi
     s_the = self.s_the
     s_phi = self.s_phi
     r11 = -c_psi*s_phi
     r12 = -c_psi*c_phi
     r13 =  0.0
     r21 =  c_the*c_phi - s_the*s_psi*s_phi
     r22 = -c_the*s_phi - s_the*s_psi*c_phi
     r23 =  0.0
     r31 =  s_the*c_phi + c_the*s_psi*s_phi
     r32 = -s_the*s_phi + c_the*s_psi*c_phi
     r33 =  0.0
     rm = matrix.sqr((r11,r12,r13, r21,r22,r23, r31,r32,r33))
     return rm
