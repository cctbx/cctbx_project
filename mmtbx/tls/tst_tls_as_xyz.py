from __future__ import division
from scitbx import matrix
from mmtbx.tls import tls_as_xyz
from libtbx import group_args
from libtbx.test_utils import approx_equal
import mmtbx.tls.tools
import iotbx.pdb
from scitbx.array_family import flex
import sys

def print_step(m, log):
  print >> log, "-"*80
  tls_as_xyz.print_step(m, log)

def extract(s):
  def fetch_matrix_1_sym(x,y,z):
    return matrix.sym(sym_mat3=[x[0], y[1], z[2], x[1], x[2], y[2]])
  def fetch_matrix_2_sym(x,y,z):
    return matrix.sym(sym_mat3=[x[3], y[4], z[5], x[4], x[5], y[5]])
  def fetch_matrix_3_sym(x,y,z):
    return matrix.sym(sym_mat3=[x[6], y[7], z[8], x[7], x[8], y[8]])
  def fetch_matrix_3_sqr(x,y,z):
    return matrix.sqr([x[6], x[7], x[8], y[6], y[7], y[8], z[6], z[7], z[8]])
  def get_next_3_lines(i, lines):
    x = lines[i+1].split()
    y = lines[i+2].split()
    z = lines[i+3].split()
    return x,y,z
  def get_next_3_lines_all_float(i, lines):
    x = [float(j) for j in lines[i+1].split()]
    y = [float(j) for j in lines[i+2].split()]
    z = [float(j) for j in lines[i+3].split()]
    return x,y,z
  read            = False
  read_T_L_S      = False
  read_lx_ly_lz   = False
  read_TL_LL_SL   = False
  read_wL         = False
  read_CLW_CLS_CL = False
  read_VL         = False
  read_VM         = False
  read_vx_vy_vz   = False
  read_VV         = False
  lines = s.splitlines()
  for i, l in enumerate(lines):
    if(l.count("************* INFORMATION FOR COMPARISON **********")):
      read=True
    if(read):
      ### set conditions
      #
      if(l.startswith("***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***")):
        read_T_L_S = True
      #
      if(l.startswith("***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)")):
        read_lx_ly_lz = True
      #
      if(l.startswith("***  T[L] L[L] S[L] *** total TLS matrices in the L base ***")):
        read_TL_LL_SL = True
      #
      if(l.startswith("***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base")):
        read_wL = True
      #
      if(l.startswith("***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***")):
        read_CLW_CLS_CL = True
      #
      if(l.startswith("***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base")):
        read_VL = True
      #
      if(l.startswith("***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base")):
        read_VM = True
      #
      if(l.startswith("***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)")):
        read_vx_vy_vz = True
      #
      if(l.startswith("***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base")):
        read_VV = True
      #
      ### extract
      #
      if(read_T_L_S):
        read_T_L_S = False
        x,y,z = get_next_3_lines_all_float(i, lines)
        T_M = fetch_matrix_1_sym(x,y,z)
        L_M = fetch_matrix_2_sym(x,y,z)
        S_M = fetch_matrix_3_sqr(x,y,z)
      #
      if(read_lx_ly_lz):
        read_lx_ly_lz = False
        x,y,z = get_next_3_lines(i, lines)
        l_x = matrix.col((float(x[1]), float(x[2]), float(x[3])))
        l_y = matrix.col((float(y[1]), float(y[2]), float(y[3])))
        l_z = matrix.col((float(z[1]), float(z[2]), float(z[3])))
      #
      if(read_TL_LL_SL):
        read_TL_LL_SL = False
        x,y,z = get_next_3_lines_all_float(i, lines)
        T_L = fetch_matrix_1_sym(x,y,z)
        L_L = fetch_matrix_2_sym(x,y,z)
        S_L = fetch_matrix_3_sqr(x,y,z)
      #
      if(read_wL):
        read_wL = False
        x,y,z = get_next_3_lines(i, lines)
        w_lx = matrix.col((float(x[7]), float(x[8]), float(x[9])))
        w_ly = matrix.col((float(y[7]), float(y[8]), float(y[9])))
        w_lz = matrix.col((float(z[7]), float(z[8]), float(z[9])))
      #
      if(read_CLW_CLS_CL):
        read_CLW_CLS_CL = False
        x,y,z = get_next_3_lines_all_float(i, lines)
        C_LW = fetch_matrix_1_sym(x,y,z)
        C_LS = fetch_matrix_2_sym(x,y,z)
        C_L  = fetch_matrix_3_sym(x,y,z)
      #
      if(read_VL):
        read_VL = False
        x,y,z = get_next_3_lines_all_float(i, lines)
        V_L = fetch_matrix_1_sym(x,y,z)
      #
      if(read_VM):
        read_VM = False
        x,y,z = get_next_3_lines_all_float(i, lines)
        V_M = fetch_matrix_1_sym(x,y,z)
      #
      if(read_vx_vy_vz):
        read_vx_vy_vz = False
        x,y,z = get_next_3_lines(i, lines)
        v_x = matrix.col((float(x[1]), float(x[2]), float(x[3])))
        v_y = matrix.col((float(y[1]), float(y[2]), float(y[3])))
        v_z = matrix.col((float(z[1]), float(z[2]), float(z[3])))
      #
      if(read_VV):
        read_VV = False
        x,y,z = get_next_3_lines_all_float(i, lines)
        V_V = fetch_matrix_1_sym(x,y,z)
      #
  return group_args(
    T_M=T_M,
    L_M=L_M,
    S_M=S_M,
    l_x=l_x,
    l_y=l_y,
    l_z=l_z,
    T_L=T_L,
    L_L=L_L,
    S_L=S_L,
    w_lx=w_lx,
    w_ly=w_ly,
    w_lz=w_lz,
    C_LW=C_LW,
    C_LS=C_LS,
    C_L =C_L ,
    V_L=V_L,
    V_M=V_M,
    v_x=v_x,
    v_y=v_y,
    v_z=v_z,
    V_V=V_V)

def compare(e, r, eps = 1.e-5):
  # e - extract object, r - result object
  assert approx_equal(e.T_M , r.T       , eps)
  assert approx_equal(e.L_M , r.L       , eps)
  assert approx_equal(e.S_M , r.S       , eps)
  assert approx_equal(e.l_x , r.b_o.l_x , eps)
  assert approx_equal(e.l_y , r.b_o.l_y , eps)
  assert approx_equal(e.l_z , r.b_o.l_z , eps)
  assert approx_equal(e.T_L , r.b_o.T_L , eps)
  assert approx_equal(e.L_L , r.b_o.L_L , eps)
  assert approx_equal(e.S_L , r.b_o.S_L , eps)
  assert approx_equal([e.w_lx[1],e.w_lx[2]], [r.c_o.w_lx[1],r.c_o.w_lx[2]], eps)
  assert approx_equal([e.w_ly[0],e.w_ly[2]], [r.c_o.w_ly[0],r.c_o.w_ly[2]], eps)
  assert approx_equal([e.w_lz[0],e.w_lz[1]], [r.c_o.w_lz[0],r.c_o.w_lz[1]], eps)
  assert approx_equal(e.C_LW, r.C_LW    , eps)
  assert approx_equal(e.C_LS, r.g_o.C_LS, eps)
  assert approx_equal(e.C_L , r.g_o.C_L , eps)
  assert approx_equal(e.V_L , r.g_o.V_L , eps)
  assert approx_equal(e.V_M , r.h_o.V_M , eps)
  assert approx_equal(e.v_x , r.h_o.v_x , eps)
  assert approx_equal(e.v_y , r.h_o.v_y , eps)
  assert approx_equal(e.v_z , r.h_o.v_z , eps)
  assert approx_equal(e.V_V , r.h_o.V_V , eps)


getTLS3_test001 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  dx  dy  dz      *** rms : Libration around lx,ly,lz
 0.1000000 0.2000000 0.3000000

***  dx2 dy2 dz2     *** rms2: Libration around lx,ly,lz
 0.0100000 0.0400000 0.0900000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005
"""

getTLS3_test004 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   2.00000   3.00000
   3.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   0.26726   0.53452   0.80178
V2=   0.92212   0.09969  -0.37383
V3=  -0.27975   0.83925  -0.46625

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.18685097 -0.13412423  0.03046584
 -0.13412423  0.45522982 -0.25211182
  0.03046584 -0.25211182  0.16791928


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.18685097 -0.13412423  0.03046584    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.13412423  0.45522982 -0.25211182    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.03046584 -0.25211182  0.16791928    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.18685097 -0.13412423  0.03046584    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.13412423  0.45522982 -0.25211182    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.03046584 -0.25211182  0.16791928    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.18685097 -0.13412423  0.03046584    0.18685097 -0.13412423  0.03046584
 -0.13412423  0.45522982 -0.25211182   -0.13412423  0.45522982 -0.25211182
  0.03046584 -0.25211182  0.16791928    0.03046584 -0.25211182  0.16791928

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.18685097 -0.13412423  0.03046584    0.18685097 -0.13412423  0.03046584
 -0.13412423  0.45522982 -0.25211182   -0.13412423  0.45522982 -0.25211182
  0.03046584 -0.25211182  0.16791928    0.03046584 -0.25211182  0.16791928

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.26726   0.53452   0.80178
Vy=   0.92212   0.09969  -0.37383
Vz=  -0.27975   0.83925  -0.46625

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.00999999  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000004 -0.00000001
  0.00000000  0.00000000  0.64000005   -0.00000001  0.00000000  0.64000005
"""

getTLS3_test011 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

principal Libration axes (orthonormal L base)
L1=   0.70711   0.70711   0.00000
L2=   0.00000   0.00000   1.00000
L3=   0.70711  -0.70711   0.00000

principal Vibration axes (orthonormal V base)
V1=   0.70711   0.70711   0.00000
V2=   0.00000   0.00000   1.00000
V3=   0.70711  -0.70711   0.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.05000000 -0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.04000000  0.05000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.32500002 -0.31500000  0.00000000
 -0.31500000  0.32500002  0.00000000
  0.00000000  0.00000000  0.16000001


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.32500002 -0.31500000  0.00000000    0.05000000 -0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.31500000  0.32500002  0.00000000   -0.04000000  0.05000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.70711   0.70711   0.00000
Ly=   0.00000   0.00000   1.00000
Lz=   0.70711  -0.70711   0.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.01000003  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000001  0.00000000  0.63999999    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.01000003  0.00000000  0.00000000    0.01000003  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000001  0.00000000  0.63999999    0.00000001  0.00000000  0.63999999

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.32500002 -0.31500000  0.00000000    0.32499999 -0.31499997  0.00000000
 -0.31500000  0.32500002  0.00000000   -0.31499997  0.32499999  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.16000001

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.70711   0.70711   0.00000
Vy=   0.00000   0.00000   1.00000
Vz=   0.70711  -0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000002  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.63999993
"""

getTLS3_test014 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   2.00000   3.00000
   3.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   1.00000
   1.00000   0.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.26726   0.53452   0.80178
L2=   0.92212   0.09969  -0.37383
L3=  -0.27975   0.83925  -0.46625

principal Vibration axes (orthonormal V base)
V1=   0.57735   0.57735   0.57735
V2=   0.81650  -0.40825  -0.40825
V3=   0.00000   0.70711  -0.70711

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.04177019 -0.01602484  0.00009317    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.01602484  0.06664596 -0.03242236    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00009317 -0.03242236  0.03158385    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.11000001 -0.05000000 -0.05000000
 -0.05000000  0.34999996 -0.28999996
 -0.05000000 -0.28999996  0.34999996


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.11000001 -0.05000000 -0.05000000    0.04177019 -0.01602484  0.00009317    0.00000000  0.00000000  0.00000000
 -0.05000000  0.34999996 -0.28999996   -0.01602484  0.06664596 -0.03242236    0.00000000  0.00000000  0.00000000
 -0.05000000 -0.28999996  0.34999996    0.00009317 -0.03242236  0.03158385    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.26726   0.53452   0.80178
Ly=   0.92212   0.09969  -0.37383
Lz=  -0.27975   0.83925  -0.46625

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.04857142 -0.08232684 -0.09121538    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.08232684  0.19281989  0.14534362    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.09121536  0.14534362  0.56860864    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.04857142 -0.08232684 -0.09121538    0.04857142 -0.08232684 -0.09121538
 -0.08232684  0.19281989  0.14534362   -0.08232684  0.19281989  0.14534362
 -0.09121536  0.14534362  0.56860864   -0.09121536  0.14534362  0.56860864

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.11000001 -0.05000000 -0.05000000    0.11000003 -0.04999999 -0.05000001
 -0.05000000  0.34999996 -0.28999996   -0.04999999  0.34999999 -0.28999996
 -0.05000000 -0.28999996  0.34999996   -0.05000002 -0.28999996  0.34999996

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.57735   0.57735   0.57735
Vy=   0.81650  -0.40825  -0.40825
Vz=   0.00000   0.70711  -0.70711

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000001  0.00000001  0.00000002
  0.00000000  0.16000001  0.00000000    0.00000001  0.16000003  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000003  0.00000001  0.63999981
"""

getTLS3_test021 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 parallel to k : (wkx,wky,wkz)=    4.00000   1.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

principal Libration axes (orthonormal L base)
L1=   0.70711   0.70711   0.00000
L2=   0.00000   0.00000   1.00000
L3=   0.70711  -0.70711   0.00000

principal Vibration axes (orthonormal V base)
V1=   0.70711   0.70711   0.00000
V2=   0.00000   0.00000   1.00000
V3=   0.70711  -0.70711   0.00000

TLS matrices from Libration in the L-base
  0.25000000 -0.36000001  0.08000001    0.01000000  0.00000000  0.00000000    0.00000000  0.03000000 -0.02000000
 -0.36000001  1.53000009 -0.06000000    0.00000000  0.04000000  0.00000000   -0.08000001  0.00000000 -0.04000000
  0.08000001 -0.06000000  0.08000001    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.00000000
TLS matrices from Libration in the M-base
  0.24499999  0.08499999 -0.29698485    0.05000000 -0.04000000  0.00000000    0.03500000  0.05500000 -0.23334524
  0.08499999  0.08499999 -0.21213204   -0.04000000  0.05000000  0.00000000   -0.05500000 -0.03500000  0.27577165
 -0.29698485 -0.21213204  1.53000009    0.00000000  0.00000000  0.04000000   -0.08485281 -0.02828427  0.00000000

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.32500002 -0.31500000  0.00000000
 -0.31500000  0.32500002  0.00000000
  0.00000000  0.00000000  0.16000001


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.56999999 -0.23000000 -0.29698485    0.05000000 -0.04000000  0.00000000    0.03500000  0.05500000 -0.23334524
 -0.23000000  0.41000003 -0.21213204   -0.04000000  0.05000000  0.00000000   -0.05500000 -0.03500000  0.27577165
 -0.29698485 -0.21213204  1.69000006    0.00000000  0.00000000  0.04000000   -0.08485281 -0.02828427  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.70711   0.70711   0.00000
Ly=   0.00000   0.00000   1.00000
Lz=   0.70711  -0.70711   0.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.25999996 -0.35999998  0.07999998    0.01000000  0.00000000  0.00000000    0.00000000  0.03000000 -0.02000000
 -0.35999998  1.69000006 -0.05999999    0.00000000  0.04000000  0.00000000   -0.08000001  0.00000000 -0.04000000
  0.07999997 -0.05999999  0.71999997    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    4.00000   1.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    2.12132  -2.12132   2.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.70711  -2.12132   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    2.82843   2.82843   1.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.25000000 -0.36000001  0.08000001    0.00000000  0.00000000  0.00000000    0.25000000 -0.36000001  0.08000001
 -0.36000001  1.53000009 -0.06000000    0.00000000  0.00000000  0.00000000   -0.36000001  1.53000009 -0.06000000
  0.08000001 -0.06000000  0.08000001    0.00000000  0.00000000  0.00000000    0.08000001 -0.06000000  0.08000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.01000003  0.00000000  0.00000000    0.00999996  0.00000003 -0.00000003
  0.00000000  0.16000001  0.00000000    0.00000003  0.15999997  0.00000001
  0.00000001  0.00000000  0.63999999   -0.00000004  0.00000001  0.63999999

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.32500002 -0.31500000  0.00000000    0.32499993 -0.31500000  0.00000003
 -0.31500000  0.32500002  0.00000000   -0.31500000  0.32500002  0.00000002
  0.00000000  0.00000000  0.16000001    0.00000003  0.00000002  0.15999997

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.70711   0.70711   0.00000
Vy=   0.00000   0.00000   1.00000
Vz=   0.70711  -0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.00999999  0.00000003 -0.00000004
  0.00000000  0.16000001  0.00000000    0.00000003  0.15999997  0.00000001
  0.00000000  0.00000000  0.64000005   -0.00000002  0.00000001  0.63999993
"""

getTLS3_test024 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   2.00000   3.00000
   3.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 parallel to k : (wkx,wky,wkz)=    4.00000   1.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   1.00000
   1.00000   0.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.26726   0.53452   0.80178
L2=   0.92212   0.09969  -0.37383
L3=  -0.27975   0.83925  -0.46625

principal Vibration axes (orthonormal V base)
V1=   0.57735   0.57735   0.57735
V2=   0.81650  -0.40825  -0.40825
V3=   0.00000   0.70711  -0.70711

TLS matrices from Libration in the L-base
  0.25000000 -0.36000001  0.08000001    0.01000000  0.00000000  0.00000000    0.00000000  0.03000000 -0.02000000
 -0.36000001  1.53000009 -0.06000000    0.00000000  0.04000000  0.00000000   -0.08000001  0.00000000 -0.04000000
  0.08000001 -0.06000000  0.08000001    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.00000000
TLS matrices from Libration in the M-base
  1.16664422 -0.06823526 -0.70200294    0.04177019 -0.01602484  0.00009317    0.08563004 -0.07749246 -0.10029086
 -0.06823528  0.16635178  0.11748905   -0.01602484  0.06664596 -0.03242236   -0.24165380 -0.00472804  0.16796263
 -0.70200288  0.11748905  0.52700424    0.00009317 -0.03242236  0.03158385    0.17404011  0.01177818 -0.08090200

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.11000001 -0.05000000 -0.05000000
 -0.05000000  0.34999996 -0.28999996
 -0.05000000 -0.28999996  0.34999996


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  1.27664423 -0.11823526 -0.75200295    0.04177019 -0.01602484  0.00009317    0.08563004 -0.07749246 -0.10029086
 -0.11823528  0.51635176 -0.17251091   -0.01602484  0.06664596 -0.03242236   -0.24165380 -0.00472804  0.16796263
 -0.75200289 -0.17251092  0.87700421    0.00009317 -0.03242236  0.03158385    0.17404011  0.01177818 -0.08090200

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.26726   0.53452   0.80178
Ly=   0.92212   0.09969  -0.37383
Lz=  -0.27975   0.83925  -0.46625

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.29857135 -0.44232675 -0.01121536    0.01000000  0.00000000  0.00000000    0.00000000  0.03000000 -0.02000000
 -0.44232675  1.72282040  0.08534358    0.00000000  0.04000000  0.00000000   -0.08000001  0.00000000 -0.04000000
 -0.01121536  0.08534360  0.64860868    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    4.00000   1.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    1.00499   2.71714  -2.14642
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -0.82676   1.14399  -1.73429
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.99117   2.23778   2.83330

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.25000000 -0.36000001  0.08000001    0.00000000  0.00000000  0.00000000    0.25000000 -0.36000001  0.08000001
 -0.36000001  1.53000009 -0.06000000    0.00000000  0.00000000  0.00000000   -0.36000001  1.53000009 -0.06000000
  0.08000001 -0.06000000  0.08000001    0.00000000  0.00000000  0.00000000    0.08000001 -0.06000000  0.08000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.04857142 -0.08232684 -0.09121538    0.04857135 -0.08232674 -0.09121536
 -0.08232684  0.19281989  0.14534362   -0.08232674  0.19282031  0.14534359
 -0.09121536  0.14534362  0.56860864   -0.09121536  0.14534360  0.56860870

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.11000001 -0.05000000 -0.05000000    0.11000045 -0.04999994 -0.05000009
 -0.05000000  0.34999996 -0.28999996   -0.04999993  0.34999999 -0.29000005
 -0.05000000 -0.28999996  0.34999996   -0.05000011 -0.29000002  0.34999990

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.57735   0.57735   0.57735
Vy=   0.81650  -0.40825  -0.40825
Vz=   0.00000   0.70711  -0.70711

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000007  0.00000024  0.00000011
  0.00000000  0.16000001  0.00000000    0.00000024  0.16000026  0.00000005
  0.00000000  0.00000000  0.64000005    0.00000010  0.00000008  0.63999987
"""

getTLS3_test031 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.500000 -0.300000  0.700000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

modified shifts sx,sy,sz for the libration axes and the trace
 -1.366667 -0.766667  0.492593            0.056000  0.000000

principal Libration axes (orthonormal L base)
L1=   0.70711   0.70711   0.00000
L2=   0.00000   0.00000   1.00000
L3=   0.70711  -0.70711   0.00000

principal Vibration axes (orthonormal V base)
V1=   0.70711   0.70711   0.00000
V2=   0.00000   0.00000   1.00000
V3=   0.70711  -0.70711   0.00000

TLS matrices from Libration in the L-base
  0.01867778  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000   -0.01366667  0.00000000  0.00000000
  0.00000000  0.02351111  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000 -0.03066667  0.00000000
  0.00000000  0.00000000  0.02183827    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.04433334
TLS matrices from Libration in the M-base
  0.02025802 -0.00158024  0.00000000    0.05000000 -0.04000000  0.00000000    0.01533333 -0.02900000  0.00000000
 -0.00158024  0.02025802  0.00000000   -0.04000000  0.05000000  0.00000000   -0.02900000  0.01533333  0.00000000
  0.00000000  0.00000000  0.02351111    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000 -0.03066667

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.32500002 -0.31500000  0.00000000
 -0.31500000  0.32500002  0.00000000
  0.00000000  0.00000000  0.16000001


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.34525806 -0.31658024  0.00000000    0.05000000 -0.04000000  0.00000000    0.01533333 -0.02900000  0.00000000
 -0.31658024  0.34525806  0.00000000   -0.04000000  0.05000000  0.00000000   -0.02900000  0.01533333  0.00000000
  0.00000000  0.00000000  0.18351112    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000 -0.03066667

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.70711   0.70711   0.00000
Ly=   0.00000   0.00000   1.00000
Lz=   0.70711  -0.70711   0.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.02867781  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000   -0.01366667  0.00000000  0.00000000
  0.00000000  0.18351112  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000 -0.03066667  0.00000000
  0.00000001  0.00000000  0.66183829    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.04433334

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
 -1.366667 -0.766667  0.492593

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.01867778  0.00000000  0.00000000    0.01867778  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.02351111  0.00000000    0.00000000  0.02351111  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.02183827    0.00000000  0.00000000  0.02183827

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.01000003  0.00000000  0.00000000    0.01000003  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000001  0.00000000  0.63999999    0.00000001  0.00000000  0.64000005

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.32500002 -0.31500000  0.00000000    0.32500002 -0.31500000  0.00000000
 -0.31500000  0.32500002  0.00000000   -0.31500000  0.32500002  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.16000001

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.70711   0.70711   0.00000
Vy=   0.00000   0.00000   1.00000
Vz=   0.70711  -0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000003  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000001  0.00000000  0.63999999
"""

getTLS3_test034 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   2.00000   3.00000
   3.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.500000 -0.300000  0.700000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   1.00000
   1.00000   0.00000   0.00000

modified shifts sx,sy,sz for the libration axes and the trace
 -1.366667 -0.766667  0.492593            0.056000  0.000000

principal Libration axes (orthonormal L base)
L1=   0.26726   0.53452   0.80178
L2=   0.92212   0.09969  -0.37383
L3=  -0.27975   0.83925  -0.46625

principal Vibration axes (orthonormal V base)
V1=   0.57735   0.57735   0.57735
V2=   0.81650  -0.40825  -0.40825
V3=   0.00000   0.70711  -0.70711

TLS matrices from Libration in the L-base
  0.01867778  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000   -0.01366667  0.00000000  0.00000000
  0.00000000  0.02351111  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000 -0.03066667  0.00000000
  0.00000000  0.00000000  0.02183827    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.04433334
TLS matrices from Libration in the M-base
  0.02303496 -0.00029772 -0.00125391    0.04177019 -0.01602484  0.00009317   -0.02358282 -0.01518012  0.01342547
 -0.00029772  0.02095190 -0.00141684   -0.01602484  0.06664596 -0.03242236   -0.01518013  0.02701657 -0.02206211
 -0.00125391 -0.00141684  0.02004031    0.00009317 -0.03242236  0.03158385    0.01342547 -0.02206211 -0.00343375

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.11000001 -0.05000000 -0.05000000
 -0.05000000  0.34999996 -0.28999996
 -0.05000000 -0.28999996  0.34999996


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.13303497 -0.05029772 -0.05125391    0.04177019 -0.01602484  0.00009317   -0.02358282 -0.01518012  0.01342547
 -0.05029772  0.37095186 -0.29141679   -0.01602484  0.06664596 -0.03242236   -0.01518013  0.02701657 -0.02206211
 -0.05125391 -0.29141679  0.37004027    0.00009317 -0.03242236  0.03158385    0.01342547 -0.02206211 -0.00343375

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.26726   0.53452   0.80178
Ly=   0.92212   0.09969  -0.37383
Lz=  -0.27975   0.83925  -0.46625

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.06724921 -0.08232684 -0.09121537    0.01000000  0.00000000  0.00000000   -0.01366667  0.00000000  0.00000000
 -0.08232683  0.21633101  0.14534362    0.00000000  0.04000000  0.00000000    0.00000000 -0.03066667  0.00000000
 -0.09121535  0.14534362  0.59044695    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.04433334

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
 -1.366667 -0.766667  0.492593

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.01867778  0.00000000  0.00000000    0.01867778  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.02351111  0.00000000    0.00000000  0.02351111  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.02183827    0.00000000  0.00000000  0.02183827

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.04857142 -0.08232684 -0.09121538    0.04857143 -0.08232684 -0.09121537
 -0.08232684  0.19281989  0.14534362   -0.08232683  0.19281989  0.14534362
 -0.09121536  0.14534362  0.56860864   -0.09121535  0.14534362  0.56860870

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.11000001 -0.05000000 -0.05000000    0.11000003 -0.05000000 -0.05000000
 -0.05000000  0.34999996 -0.28999996   -0.05000000  0.35000002 -0.28999996
 -0.05000000 -0.28999996  0.34999996   -0.05000002 -0.28999999  0.34999996

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.57735   0.57735   0.57735
Vy=   0.81650  -0.40825  -0.40825
Vz=   0.00000   0.70711  -0.70711

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000002  0.00000000  0.00000001
  0.00000000  0.16000001  0.00000000    0.00000001  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000006 -0.00000001  0.63999981
"""

getTLS3_test041 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 parallel to k : (wkx,wky,wkz)=    4.00000   1.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.500000 -0.300000  0.700000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   0.00000
   0.00000   0.00000   1.00000

modified shifts sx,sy,sz for the libration axes and the trace
 -1.366667 -0.766667  0.492593            0.056000  0.000000

principal Libration axes (orthonormal L base)
L1=   0.70711   0.70711   0.00000
L2=   0.00000   0.00000   1.00000
L3=   0.70711  -0.70711   0.00000

principal Vibration axes (orthonormal V base)
V1=   0.70711   0.70711   0.00000
V2=   0.00000   0.00000   1.00000
V3=   0.70711  -0.70711   0.00000

TLS matrices from Libration in the L-base
  0.26867777 -0.36000001  0.08000001    0.01000000  0.00000000  0.00000000   -0.01366667  0.03000000 -0.02000000
 -0.36000001  1.55351114 -0.06000000    0.00000000  0.04000000  0.00000000   -0.08000001 -0.03066667 -0.04000000
  0.08000001 -0.06000000  0.10183828    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.04433334
TLS matrices from Libration in the M-base
  0.26525804  0.08341976 -0.29698485    0.05000000 -0.04000000  0.00000000    0.05033333  0.02600000 -0.23334524
  0.08341974  0.10525802 -0.21213204   -0.04000000  0.05000000  0.00000000   -0.08400000 -0.01966666  0.27577165
 -0.29698485 -0.21213204  1.55351114    0.00000000  0.00000000  0.04000000   -0.08485281 -0.02828427 -0.03066667

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.32500002 -0.31500000  0.00000000
 -0.31500000  0.32500002  0.00000000
  0.00000000  0.00000000  0.16000001


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.59025806 -0.23158024 -0.29698485    0.05000000 -0.04000000  0.00000000    0.05033333  0.02600000 -0.23334524
 -0.23158026  0.43025804 -0.21213204   -0.04000000  0.05000000  0.00000000   -0.08400000 -0.01966666  0.27577165
 -0.29698485 -0.21213204  1.71351111    0.00000000  0.00000000  0.04000000   -0.08485281 -0.02828427 -0.03066667

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.70711   0.70711   0.00000
Ly=   0.00000   0.00000   1.00000
Lz=   0.70711  -0.70711   0.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.27867779 -0.35999998  0.08000001    0.01000000  0.00000000  0.00000000   -0.01366667  0.03000000 -0.02000000
 -0.35999998  1.71351111 -0.05999999    0.00000000  0.04000000  0.00000000   -0.08000001 -0.03066667 -0.04000000
  0.08000002 -0.05999999  0.74183828    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.04433334

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    4.00000   1.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    2.12132  -2.12132   2.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.70711  -2.12132   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    2.82843   2.82843   1.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
 -1.366667 -0.766667  0.492593

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.25000000 -0.36000001  0.08000001    0.01867778  0.00000000  0.00000000    0.26867777 -0.36000001  0.08000001
 -0.36000001  1.53000009 -0.06000000    0.00000000  0.02351111  0.00000000   -0.36000001  1.55351114 -0.06000000
  0.08000001 -0.06000000  0.08000001    0.00000000  0.00000000  0.02183827    0.08000001 -0.06000000  0.10183828

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.01000003  0.00000000  0.00000000    0.01000002  0.00000003  0.00000001
  0.00000000  0.16000001  0.00000000    0.00000003  0.15999997  0.00000001
  0.00000001  0.00000000  0.63999999    0.00000001  0.00000001  0.63999999

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.32500002 -0.31500000  0.00000000    0.32499999 -0.31499997  0.00000003
 -0.31500000  0.32500002  0.00000000   -0.31499997  0.32499999  0.00000002
  0.00000000  0.00000000  0.16000001    0.00000003  0.00000002  0.15999997

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.70711   0.70711   0.00000
Vy=   0.00000   0.00000   1.00000
Vz=   0.70711  -0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000002  0.00000003  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000003  0.15999997  0.00000001
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000001  0.63999993
"""

getTLS3_test044 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   2.00000   3.00000
   3.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 parallel to k : (wkx,wky,wkz)=    4.00000   1.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.500000 -0.300000  0.700000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.1000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0100000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000   1.00000
   1.00000   0.00000   0.00000

modified shifts sx,sy,sz for the libration axes and the trace
 -1.366667 -0.766667  0.492593            0.056000  0.000000

principal Libration axes (orthonormal L base)
L1=   0.26726   0.53452   0.80178
L2=   0.92212   0.09969  -0.37383
L3=  -0.27975   0.83925  -0.46625

principal Vibration axes (orthonormal V base)
V1=   0.57735   0.57735   0.57735
V2=   0.81650  -0.40825  -0.40825
V3=   0.00000   0.70711  -0.70711

TLS matrices from Libration in the L-base
  0.26867777 -0.36000001  0.08000001    0.01000000  0.00000000  0.00000000   -0.01366667  0.03000000 -0.02000000
 -0.36000001  1.55351114 -0.06000000    0.00000000  0.04000000  0.00000000   -0.08000001 -0.03066667 -0.04000000
  0.08000001 -0.06000000  0.10183828    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.04433334
TLS matrices from Libration in the M-base
  1.18967915 -0.06853300 -0.70325685    0.04177019 -0.01602484  0.00009317    0.06204721 -0.09267259 -0.08686541
 -0.06853301  0.18730368  0.11607221   -0.01602484  0.06664596 -0.03242236   -0.25683391  0.02228853  0.14590052
 -0.70325685  0.11607219  0.54704458    0.00009317 -0.03242236  0.03158385    0.18746558 -0.01028393 -0.08433574

V matrix from Vibration in the V-base
  0.01000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.11000001 -0.05000000 -0.05000000
 -0.05000000  0.34999996 -0.28999996
 -0.05000000 -0.28999996  0.34999996


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  1.29967916 -0.11853299 -0.75325686    0.04177019 -0.01602484  0.00009317    0.06204721 -0.09267259 -0.08686541
 -0.11853301  0.53730363 -0.17392775   -0.01602484  0.06664596 -0.03242236   -0.25683391  0.02228853  0.14590052
 -0.75325686 -0.17392777  0.89704454    0.00009317 -0.03242236  0.03158385    0.18746558 -0.01028393 -0.08433574

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.26726   0.53452   0.80178
Ly=   0.92212   0.09969  -0.37383
Lz=  -0.27975   0.83925  -0.46625

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.31724912 -0.44232687 -0.01121537    0.01000000  0.00000000  0.00000000   -0.01366667  0.03000000 -0.02000000
 -0.44232678  1.74633145  0.08534361    0.00000000  0.04000000  0.00000000   -0.08000001 -0.03066667 -0.04000000
 -0.01121539  0.08534359  0.67044687    0.00000000  0.00000000  0.09000000    0.09000000 -0.36000001  0.04433334

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   2.00000   3.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   0.00000   2.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    4.00000   1.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    1.00499   2.71714  -2.14642
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -0.82676   1.14399  -1.73429
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.99117   2.23778   2.83330

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
 -1.366667 -0.766667  0.492593

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.25000000 -0.36000001  0.08000001    0.01867778  0.00000000  0.00000000    0.26867777 -0.36000001  0.08000001
 -0.36000001  1.53000009 -0.06000000    0.00000000  0.02351111  0.00000000   -0.36000001  1.55351114 -0.06000000
  0.08000001 -0.06000000  0.08000001    0.00000000  0.00000000  0.02183827    0.08000001 -0.06000000  0.10183828

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.04857142 -0.08232684 -0.09121538    0.04857135 -0.08232686 -0.09121538
 -0.08232684  0.19281989  0.14534362   -0.08232677  0.19282031  0.14534362
 -0.09121536  0.14534362  0.56860864   -0.09121539  0.14534360  0.56860858

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.11000001 -0.05000000 -0.05000000    0.11000040 -0.04999991 -0.05000012
 -0.05000000  0.34999996 -0.28999996   -0.04999998  0.34999990 -0.28999999
 -0.05000000 -0.28999996  0.34999996   -0.05000018 -0.28999999  0.34999993

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.57735   0.57735   0.57735
Vy=   0.81650  -0.40825  -0.40825
Vz=   0.00000   0.70711  -0.70711

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.01000000  0.00000000  0.00000000    0.01000003  0.00000015  0.00000007
  0.00000000  0.16000001  0.00000000    0.00000025  0.16000028  0.00000013
  0.00000000  0.00000000  0.64000005    0.00000006  0.00000013  0.63999975
"""

# Ensemble generation tests ####################################################

pdb_str = """\n
REMARK this is PDB formatted CH3-HN-CO-CH3-mod.xyz file
CRYST1    6.000    6.000    3.000  90.00  90.00  90.00    P 1
ATOM      1  N    ALA    1       1.500   0.000   0.000  1.00  0.91      E
ATOM      1  CA   ALA    1       0.000   0.000   0.000  1.00  0.91      E
ATOM      1  CB   ALA    1       0.000   1.500   0.000  1.00  0.91      E
ATOM      1  C    ALA    1       0.000   0.000   1.500  1.00  0.91      E
ATOM      1  O    ALA    1      -1.000  -1.000   2.000  1.00  0.91      E
end
"""

getTLS_test111="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test112="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   1.00000  -1.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   1.00000  -1.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   1.00000  -1.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.04000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test113="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test114 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   1.00000  -1.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001 -0.04000000 -0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   1.00000  -1.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   1.00000  -1.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.04000000
  0.00000000  0.04000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.04000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test121 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.2000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0400000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0400000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.2000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test122 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.2000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0400000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    1.00000   0.00000  -1.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.00000000  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.00000000  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.00000000  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.00000000  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0400000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.2000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.00000   0.00000  -1.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.00000   0.00000  -1.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.04000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.04000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test123 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.2000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0400000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  2.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.08000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.08000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.08000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.08000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0400000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.2000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  2.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test124 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.2000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0400000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    1.00000   0.00000  -1.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  2.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.08000001  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.08000001  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.08000001  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.04000000  0.00000000    0.04000000  0.08000001  0.04000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0400000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.2000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.00000   0.00000  -1.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.00000   0.00000  -1.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  2.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.04000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.04000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.04000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test131 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test132 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    1.00000  -1.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.00000000
TLS matrices from Libration in the M-base
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.00000  -1.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.00000  -1.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.04000000  0.04000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.04000000  0.04000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test133 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.08000001
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.08000001

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.08000001

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.08000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.16000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test134 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    1.00000  -1.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.08000001
TLS matrices from Libration in the M-base
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.08000001

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.08000001

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000   -0.04000000 -0.04000000  0.08000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.00000  -1.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.00000  -1.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.04000000  0.04000000  0.00000000
  0.04000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000    0.04000000  0.04000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.16000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test141="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.04000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.04000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005
"""

getTLS_test142="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   0.66667   0.66667  -0.33333
V2=  -0.70711   0.70711   0.00000
V3=   0.23570   0.23570   0.94281

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.13333334 -0.02666666  0.13333334
 -0.02666666  0.13333334  0.13333334
  0.13333334  0.13333334  0.57333338


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.13333334 -0.02666666  0.13333334    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.02666666  0.13333334  0.13333334    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.13333334  0.13333334  0.57333338    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.13333334 -0.02666666  0.13333334    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.02666666  0.13333334  0.13333334    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.13333334  0.13333334  0.57333338    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.13333334 -0.02666666  0.13333334    0.13333334 -0.02666666  0.13333334
 -0.02666666  0.13333334  0.13333334   -0.02666666  0.13333334  0.13333334
  0.13333334  0.13333334  0.57333338    0.13333334  0.13333334  0.57333338

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.13333334 -0.02666666  0.13333334    0.13333334 -0.02666666  0.13333334
 -0.02666666  0.13333334  0.13333334   -0.02666666  0.13333334  0.13333334
  0.13333334  0.13333334  0.57333338    0.13333334  0.13333334  0.57333338

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.66667   0.66667  -0.33333
Vy=  -0.70711   0.70711   0.00000
Vz=   0.23570   0.23570   0.94281

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.04000001  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000000  0.00000000
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005
"""

getTLS_test143="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
  -1.00000   1.00000   0.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=  -0.70711   0.70711   0.00000
V2=   0.66667   0.66667  -0.33333
V3=  -0.23570  -0.23570  -0.94281

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.12666668  0.08666668  0.10666668
  0.08666668  0.12666668  0.10666668
  0.10666668  0.10666668  0.58666670


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.12666668  0.08666668  0.10666668    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.08666668  0.12666668  0.10666668    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.10666668  0.10666668  0.58666670    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.12666668  0.08666668  0.10666668    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.08666668  0.12666668  0.10666668    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.10666668  0.10666668  0.58666670    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.12666668  0.08666668  0.10666668    0.12666668  0.08666668  0.10666668
  0.08666668  0.12666668  0.10666668    0.08666668  0.12666668  0.10666668
  0.10666668  0.10666668  0.58666670    0.10666668  0.10666668  0.58666670

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.12666668  0.08666668  0.10666668    0.12666668  0.08666668  0.10666668
  0.08666668  0.12666668  0.10666668    0.08666668  0.12666668  0.10666668
  0.10666668  0.10666668  0.58666670    0.10666668  0.10666668  0.58666670

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=  -0.70711   0.70711   0.00000
Vy=   0.66667   0.66667  -0.33333
Vz=  -0.23570  -0.23570  -0.94281

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.03999999  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000004 -0.00000001
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000000  0.64000005
"""

getTLS_test144="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the L-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.40000001 -0.24000001  0.00000000
 -0.24000001  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.40000001 -0.24000001  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.24000001  0.40000001  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.40000001 -0.24000001  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.24000001  0.40000001  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.40000001 -0.24000001  0.00000000    0.40000001 -0.24000001  0.00000000
 -0.24000001  0.40000001  0.00000000   -0.24000001  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.04000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.40000001 -0.24000001  0.00000000    0.40000001 -0.24000001  0.00000000
 -0.24000001  0.40000001  0.00000000   -0.24000001  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.04000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000000 -0.00000001
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000001  0.63999999
"""

getTLS_test161 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   0.23570   0.23570   0.94281
  -0.70711   0.70711   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.23570   0.23570   0.94281
L2=  -0.70711   0.70711   0.00000
L3=  -0.66667  -0.66667   0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.01777782  0.01777782 -0.00888882    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.01777782  0.01777782 -0.00888882    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888882 -0.00888882  0.00444436    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.01777782  0.01777782 -0.00888882    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.01777782  0.01777782 -0.00888882    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888882 -0.00888882  0.00444436    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.23570   0.23570   0.94281
Ly=  -0.70711   0.70711   0.00000
Lz=  -0.66667  -0.66667   0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test162 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   0.23570   0.23570   0.94281
  -0.70711   0.70711   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   1.50000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.23570   0.23570   0.94281
L2=  -0.70711   0.70711   0.00000
L3=  -0.66667  -0.66667   0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000021  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000 -0.05656862  0.00000000
TLS matrices from Libration in the M-base
  0.04000010 -0.04000010  0.00000000    0.01777782  0.01777782 -0.00888882   -0.02666673  0.02666673  0.00000000
 -0.04000010  0.04000010  0.00000000    0.01777782  0.01777782 -0.00888882   -0.02666673  0.02666673  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888882 -0.00888882  0.00444436    0.01333322 -0.01333322  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000010 -0.04000010  0.00000000    0.01777782  0.01777782 -0.00888882   -0.02666673  0.02666673  0.00000000
 -0.04000010  0.04000010  0.00000000    0.01777782  0.01777782 -0.00888882   -0.02666673  0.02666673  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888882 -0.00888882  0.00444436    0.01333322 -0.01333322  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.23570   0.23570   0.94281
Ly=  -0.70711   0.70711   0.00000
Lz=  -0.66667  -0.66667   0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000020  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.04000000    0.00000000 -0.05656862  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.41422   0.00000   0.50000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000021  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.08000021  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test163 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   0.23570   0.23570   0.94281
  -0.70711   0.70711   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.23570   0.23570   0.94281
L2=  -0.70711   0.70711   0.00000
L3=  -0.66667  -0.66667   0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.08000001
TLS matrices from Libration in the M-base
  0.07111128  0.07111128 -0.03555527    0.01777782  0.01777782 -0.00888882    0.03555564  0.03555564 -0.01777763
  0.07111128  0.07111128 -0.03555527    0.01777782  0.01777782 -0.00888882    0.03555564  0.03555564 -0.01777763
 -0.03555526 -0.03555526  0.01777744   -0.00888882 -0.00888882  0.00444436   -0.01777763 -0.01777763  0.00888872

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.07111128  0.07111128 -0.03555527    0.01777782  0.01777782 -0.00888882    0.03555564  0.03555564 -0.01777763
  0.07111128  0.07111128 -0.03555527    0.01777782  0.01777782 -0.00888882    0.03555564  0.03555564 -0.01777763
 -0.03555526 -0.03555526  0.01777744   -0.00888882 -0.00888882  0.00444436   -0.01777763 -0.01777763  0.00888872

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.23570   0.23570   0.94281
Ly=  -0.70711   0.70711   0.00000
Lz=  -0.66667  -0.66667   0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000001  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.08000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.16000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test164 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.2000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0400000
vectors defining the principal Libration axes
   0.23570   0.23570   0.94281
  -0.70711   0.70711   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   1.50000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.23570   0.23570   0.94281
L2=  -0.70711   0.70711   0.00000
L3=  -0.66667  -0.66667   0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000021  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000 -0.05656862  0.08000001
TLS matrices from Libration in the M-base
  0.11111139  0.03111119 -0.03555527    0.01777782  0.01777782 -0.00888882    0.00888891  0.06222238 -0.01777763
  0.03111119  0.11111139 -0.03555527    0.01777782  0.01777782 -0.00888882    0.00888891  0.06222238 -0.01777763
 -0.03555526 -0.03555526  0.01777744   -0.00888882 -0.00888882  0.00444436   -0.00444441 -0.03111086  0.00888872

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.11111139  0.03111119 -0.03555527    0.01777782  0.01777782 -0.00888882    0.00888891  0.06222238 -0.01777763
  0.03111119  0.11111139 -0.03555527    0.01777782  0.01777782 -0.00888882    0.00888891  0.06222238 -0.01777763
 -0.03555526 -0.03555526  0.01777744   -0.00888882 -0.00888882  0.00444436   -0.00444441 -0.03111086  0.00888872

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.23570   0.23570   0.94281
Ly=  -0.70711   0.70711   0.00000
Lz=  -0.66667  -0.66667   0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000 -0.00000001    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000020  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.04000000    0.00000000 -0.05656862  0.08000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0400000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.2000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    1.41422   0.00000   0.50000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000021  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.08000021  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.16000001    0.00000000  0.00000000  0.16000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000
"""

getTLS_test171 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.01777778  0.01777778 -0.00888889    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.01777778  0.01777778 -0.00888889    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888889 -0.00888889  0.00444445    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.01777778  0.01777778 -0.00888889    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.01777778  0.01777778 -0.00888889    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888889 -0.00888889  0.00444445    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test172 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.05656854  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.04000000 -0.04000000  0.00000000    0.01777778  0.01777778 -0.00888889   -0.02666667  0.02666667  0.00000000
 -0.04000000  0.04000000  0.00000000    0.01777778  0.01777778 -0.00888889   -0.02666667  0.02666667  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888889 -0.00888889  0.00444445    0.01333333 -0.01333333  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.04000000 -0.04000000  0.00000000    0.01777778  0.01777778 -0.00888889   -0.02666667  0.02666667  0.00000000
 -0.04000000  0.04000000  0.00000000    0.01777778  0.01777778 -0.00888889   -0.02666667  0.02666667  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00888889 -0.00888889  0.00444445    0.01333333 -0.01333333  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.00000000  0.05656854  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.08000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test173 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.07111112  0.07111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.03555556  0.03555556 -0.01777778
  0.07111112  0.07111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.03555556  0.03555556 -0.01777778
 -0.03555556 -0.03555556  0.01777778   -0.00888889 -0.00888889  0.00444445   -0.01777778 -0.01777778  0.00888889

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.07111112  0.07111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.03555556  0.03555556 -0.01777778
  0.07111112  0.07111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.03555556  0.03555556 -0.01777778
 -0.03555556 -0.03555556  0.01777778   -0.00888889 -0.00888889  0.00444445   -0.01777778 -0.01777778  0.00888889

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.16000004  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000003  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000001 -0.00000001
  0.00000000  0.00000000  0.00000000   -0.00000001 -0.00000001  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000001 -0.00000001
  0.00000000  0.00000000  0.00000000   -0.00000001 -0.00000001  0.00000000
"""

getTLS_test174 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.05656854  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.11111112  0.03111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
  0.03111112  0.11111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
 -0.03555556 -0.03555556  0.01777778   -0.00888889 -0.00888889  0.00444445   -0.00444445 -0.03111111  0.00888889

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.11111112  0.03111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
  0.03111112  0.11111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
 -0.03555556 -0.03555556  0.01777778   -0.00888889 -0.00888889  0.00444445   -0.00444445 -0.03111111  0.00888889

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.16000004  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.05656854  0.00000000
  0.00000000  0.07999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.08000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000003  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000002  0.00000001 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000002  0.00000001 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
"""

getTLS_test174new="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.2000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0400000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.16000001  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.05656854  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.11111112  0.03111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
  0.03111112  0.11111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
 -0.03555556 -0.03555556  0.01777778   -0.00888889 -0.00888889  0.00444445   -0.00444445 -0.03111111  0.00888889

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.11111112  0.03111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
  0.03111112  0.11111112 -0.03555556    0.01777778  0.01777778 -0.00888889    0.00888889  0.06222223 -0.01777778
 -0.03555556 -0.03555556  0.01777778   -0.00888889 -0.00888889  0.00444445   -0.00444445 -0.03111111  0.00888889

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.16000004  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000    0.08000001  0.05656854  0.00000000
  0.00000000  0.07999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0400000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.2000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000    0.16000001  0.00000000  0.00000000
  0.00000000  0.08000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.08000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000003  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000002  0.00000001 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000002  0.00000001 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
"""

getTLS_test175="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.3000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0900000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.36000001  0.00000000  0.00000000    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.25000000  0.07000002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.07000002  0.25000000 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04000000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.25000000  0.07000002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.07000002  0.25000000 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04000000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.36000004  0.00000000  0.00000000    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.17999998  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0900000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.3000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.17999999  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000003  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000002 -0.00000001
  0.00000000  0.00000000  0.00000000    0.00000002  0.00000000 -0.00000001
  0.00000000  0.00000000  0.00000000   -0.00000001  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00000001  0.00000003  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000 -0.00000001
"""

getTLS_test176="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0200000 0.0400000 0.0800000
tx2,ty2,tz2= 0.0004000 0.0016000 0.0064000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00040000  0.00000000  0.00000000
  0.00000000  0.00160000  0.00000000
  0.00000000  0.00000000  0.00640000
V matrix from Vibration in the M-base
  0.00400000 -0.00240000  0.00000000
 -0.00240000  0.00400000  0.00000000
  0.00000000  0.00000000  0.00040000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00400000 -0.00240000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
 -0.00240000  0.00400000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00040000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00146667  0.00000000  0.00037712    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00640000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00037712  0.00000000  0.00053333    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00146667  0.00000000  0.00037712    0.00146667  0.00000000  0.00037712
  0.00000000  0.00640000  0.00000000    0.00000000  0.00640000  0.00000000
  0.00037712  0.00000000  0.00053333    0.00037712  0.00000000  0.00053333

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00400000 -0.00240000  0.00000000    0.00400000 -0.00240000  0.00000000
 -0.00240000  0.00400000  0.00000000   -0.00240000  0.00400000  0.00000000
  0.00000000  0.00000000  0.00040000    0.00000000  0.00000000  0.00040000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00040000  0.00000000  0.00000000    0.00040000  0.00000000  0.00000000
  0.00000000  0.00160000  0.00000000    0.00000000  0.00160000  0.00000000
  0.00000000  0.00000000  0.00640000    0.00000000  0.00000000  0.00640000
"""

getTLS_test181 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.3000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0900000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0200000 0.0400000 0.0800000
tx2,ty2,tz2= 0.0004000 0.0016000 0.0064000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.36000001  0.00000000  0.00000000    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.25000000  0.07000002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.07000002  0.25000000 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04000000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

V matrix from Vibration in the V-base
  0.00040000  0.00000000  0.00000000
  0.00000000  0.00160000  0.00000000
  0.00000000  0.00000000  0.00640000
V matrix from Vibration in the M-base
  0.00400000 -0.00240000  0.00000000
 -0.00240000  0.00400000  0.00000000
  0.00000000  0.00000000  0.00040000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.25400001  0.06760002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.06760002  0.25400001 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04040000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.36146674 -0.00000001  0.00037712    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.18639998  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00037713  0.00000000  0.00053333    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0900000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.3000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.17999999  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00146667  0.00000000  0.00037712    0.00146672 -0.00000001  0.00037712
  0.00000000  0.00640000  0.00000000    0.00000000  0.00639999  0.00000000
  0.00037712  0.00000000  0.00053333    0.00037713  0.00000000  0.00053333

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00400000 -0.00240000  0.00000000    0.00400002 -0.00239997 -0.00000001
 -0.00240000  0.00400000  0.00000000   -0.00239997  0.00400002 -0.00000002
  0.00000000  0.00000000  0.00040000   -0.00000001 -0.00000001  0.00040001

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00040000  0.00000000  0.00000000    0.00040001 -0.00000001  0.00000000
  0.00000000  0.00160000  0.00000000   -0.00000002  0.00160005  0.00000000
  0.00000000  0.00000000  0.00640000    0.00000000  0.00000000  0.00639999
"""

getTLS_test181new="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.3000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0900000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0200000 0.0400000 0.0800000
tx2,ty2,tz2= 0.0004000 0.0016000 0.0064000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.36000001  0.00000000  0.00000000    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.25000000  0.07000002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.07000002  0.25000000 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04000000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

V matrix from Vibration in the V-base
  0.00040000  0.00000000  0.00000000
  0.00000000  0.00160000  0.00000000
  0.00000000  0.00000000  0.00640000
V matrix from Vibration in the M-base
  0.00400000 -0.00240000  0.00000000
 -0.00240000  0.00400000  0.00000000
  0.00000000  0.00000000  0.00040000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.25400001  0.06760002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.06760002  0.25400001 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04040000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.36146674 -0.00000001  0.00037712    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.18639998  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00037713  0.00000000  0.00053333    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0900000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.3000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.17999999  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00146667  0.00000000  0.00037712    0.00146672 -0.00000001  0.00037712
  0.00000000  0.00640000  0.00000000    0.00000000  0.00639999  0.00000000
  0.00037712  0.00000000  0.00053333    0.00037713  0.00000000  0.00053333

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00400000 -0.00240000  0.00000000    0.00400002 -0.00239997 -0.00000001
 -0.00240000  0.00400000  0.00000000   -0.00239997  0.00400002 -0.00000002
  0.00000000  0.00000000  0.00040000   -0.00000001 -0.00000001  0.00040001

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00040000  0.00000000  0.00000000    0.00040001 -0.00000001  0.00000000
  0.00000000  0.00160000  0.00000000   -0.00000002  0.00160005  0.00000000
  0.00000000  0.00000000  0.00640000    0.00000000  0.00000000  0.00639999
"""

getTLS_test182 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.0300000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0009000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.00360000  0.00000000  0.00000000    0.00090000  0.00000000  0.00000000    0.00180000  0.00127279  0.00000000
  0.00000000  0.00180000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00250000  0.00070000 -0.00080000    0.00040000  0.00040000 -0.00020000    0.00020000  0.00140000 -0.00040000
  0.00070000  0.00250000 -0.00080000    0.00040000  0.00040000 -0.00020000    0.00020000  0.00140000 -0.00040000
 -0.00080000 -0.00080000  0.00040000   -0.00020000 -0.00020000  0.00010000   -0.00010000 -0.00070000  0.00020000

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.40000001 -0.24000001  0.00000000
 -0.24000001  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.40250000 -0.23930001 -0.00080000    0.00040000  0.00040000 -0.00020000    0.00020000  0.00140000 -0.00040000
 -0.23930001  0.40250000 -0.00080000    0.00040000  0.00040000 -0.00020000    0.00020000  0.00140000 -0.00040000
 -0.00080000 -0.00080000  0.04040000   -0.00020000 -0.00020000  0.00010000   -0.00010000 -0.00070000  0.00020000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.15026666 -0.00000001  0.03771236    0.00090000  0.00000000  0.00000000    0.00180000  0.00127279  0.00000000
 -0.00000003  0.64179999 -0.00000001    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.03771236  0.00000000  0.05333333    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0009000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.0300000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00360000  0.00000000  0.00000000    0.00360000  0.00000000  0.00000000
  0.00000000  0.00180000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00180000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.14666668  0.00000000  0.03771236    0.14666666 -0.00000001  0.03771236
 -0.00000001  0.63999999  0.00000000   -0.00000003  0.63999999 -0.00000001
  0.03771236  0.00000000  0.05333334    0.03771236  0.00000000  0.05333333

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.40000001 -0.24000001  0.00000000    0.40000001 -0.23999996  0.00000000
 -0.24000001  0.40000001  0.00000000   -0.23999999  0.39999998  0.00000000
  0.00000000  0.00000000  0.04000000    0.00000000  0.00000000  0.04000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000001 -0.00000001
  0.00000000  0.00000000  0.64000005    0.00000000 -0.00000003  0.63999993
"""

getTLS_test183 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.3000000 0.0000000 0.0000000
dx2,dy2,dz2= 0.0900000 0.0000000 0.0000000
vectors defining the principal Libration axes
   1.00000   1.00000  -0.50000
  -1.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=   0.66667   0.66667  -0.33333
L2=  -0.70711   0.70711   0.00000
L3=   0.23570   0.23570   0.94281

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.36000001  0.00000000  0.00000000    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.25000000  0.07000002 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
  0.07000002  0.25000000 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.04000000   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.40000001 -0.24000001  0.00000000
 -0.24000001  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.64999998 -0.16999999 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.16999999  0.64999998 -0.08000001    0.04000000  0.04000000 -0.02000000    0.02000001  0.14000000 -0.04000000
 -0.08000001 -0.08000001  0.08000001   -0.02000000 -0.02000000  0.01000000   -0.01000000 -0.07000000  0.02000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   0.66667   0.66667  -0.33333
Ly=  -0.70711   0.70711   0.00000
Lz=   0.23570   0.23570   0.94281

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.50666672 -0.00000002  0.03771237    0.09000000  0.00000000  0.00000000    0.18000001  0.12727922  0.00000000
  0.00000003  0.81999993  0.00000001    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.03771235 -0.00000001  0.05333333    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0900000   0.0000000   0.0000000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.3000000   0.0000000   0.0000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=   -0.50000   0.00000   1.41421
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   1.50000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  2.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000    0.36000001  0.00000000  0.00000000
  0.00000000  0.17999999  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.17999999  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.14666668  0.00000000  0.03771236    0.14666671 -0.00000002  0.03771237
 -0.00000001  0.63999999  0.00000000    0.00000003  0.63999993  0.00000001
  0.03771236  0.00000000  0.05333334    0.03771235 -0.00000001  0.05333333

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.40000001 -0.24000001  0.00000000    0.39999998 -0.23999996  0.00000000
 -0.24000001  0.40000001  0.00000000   -0.23999992  0.39999998 -0.00000001
  0.00000000  0.00000000  0.04000000   -0.00000001 -0.00000001  0.04000001

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.04000001 -0.00000002  0.00000000
  0.00000000  0.16000001  0.00000000    0.00000000  0.16000003 -0.00000003
  0.00000000  0.00000000  0.64000005    0.00000000  0.00000003  0.63999987
"""

getTLS_test190 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   1.50000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=   1.00000   0.00000   0.00000
L2=   0.00000   1.00000   0.00000
L3=   0.00000   0.00000   1.00000

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.18000001
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.18000001

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.18000001

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=   1.00000   0.00000   0.00000
Ly=   0.00000   1.00000   0.00000
Lz=   0.00000   0.00000   1.00000

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.18000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.36000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test191 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
  -0.70711   0.70711   0.00000
   0.23570   0.23570   0.94281

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=  -0.70711   0.70711   0.00000
L2=   0.23570   0.23570   0.94281
L3=   0.66667   0.66667  -0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.00000000  0.00000000  0.00000000    0.04722228  0.03722228 -0.01111102    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.03722228  0.04722228 -0.01111102    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.01111102 -0.01111102  0.04555546    0.00000000  0.00000000  0.00000000

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.00000000  0.00000000  0.00000000    0.04722228  0.03722228 -0.01111102    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.03722228  0.04722228 -0.01111102    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.01111102 -0.01111102  0.04555546    0.00000000  0.00000000  0.00000000

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=  -0.70711   0.70711   0.00000
Ly=   0.23570   0.23570   0.94281
Lz=   0.66667   0.66667  -0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test192 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
  -0.70711   0.70711   0.00000
   0.23570   0.23570   0.94281

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    1.00000   1.20000  -2.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   1.30000   0.50000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   1.50000
correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=  -0.70711   0.70711   0.00000
L2=   0.23570   0.23570   0.94281
L3=   0.66667   0.66667  -0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.18004490  0.00000000 -0.00216858    0.01000000  0.00000000  0.00000000    0.00000000  0.02133329  0.01367080
  0.00000000  0.04551093  0.02916432    0.00000000  0.04000000  0.00000000   -0.00133340  0.00000000  0.06505383
 -0.00216858  0.02916432  0.12448908    0.00000000  0.00000000  0.09000000    0.12727939  0.00000000  0.00000000
TLS matrices from Libration in the M-base
  0.15908936 -0.02300013 -0.00202204    0.04722228  0.03722228 -0.01111102   -0.05955578  0.06000006 -0.01611102
 -0.02300012  0.15500022 -0.00099978    0.03722228  0.04722228 -0.01111102   -0.03955578  0.08000006  0.00588899
 -0.00202205 -0.00099978  0.03595534   -0.01111102 -0.01111102  0.04555546    0.07177768  0.01000030 -0.02044428

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.15908936 -0.02300013 -0.00202204    0.04722228  0.03722228 -0.01111102   -0.05955578  0.06000006 -0.01611102
 -0.02300012  0.15500022 -0.00099978    0.03722228  0.04722228 -0.01111102   -0.03955578  0.08000006  0.00588899
 -0.00202205 -0.00099978  0.03595534   -0.01111102 -0.01111102  0.04555546    0.07177768  0.01000030 -0.02044428

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=  -0.70711   0.70711   0.00000
Ly=   0.23570   0.23570   0.94281
Lz=   0.66667   0.66667  -0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.18004490  0.00000000 -0.00216858    0.01000000  0.00000000  0.00000000    0.00000000  0.02133329  0.01367080
  0.00000000  0.04551093  0.02916432    0.00000000  0.04000000  0.00000000   -0.00133340  0.00000000  0.06505383
 -0.00216859  0.02916432  0.12448908    0.00000000  0.00000000  0.09000000    0.12727939  0.00000000  0.00000000

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.14142  -1.36708   2.13333
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.62635   0.54212   0.03334
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   1.41422  -0.50000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    1.00000   1.20000  -2.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   1.30000   0.50000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  0.000000  0.000000  0.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.18004490  0.00000000 -0.00216858    0.00000000  0.00000000  0.00000000    0.18004490  0.00000000 -0.00216858
  0.00000000  0.04551093  0.02916432    0.00000000  0.00000000  0.00000000    0.00000000  0.04551093  0.02916432
 -0.00216858  0.02916432  0.12448908    0.00000000  0.00000000  0.00000000   -0.00216858  0.02916432  0.12448908

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00000001  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000001 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000001 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
"""

getTLS_test193 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
  -0.70711   0.70711   0.00000
   0.23570   0.23570   0.94281

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 parallel to j : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   0.00000
correlation shifts sx,sy,sz for the libration axes
  1.000000 -1.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=  -0.70711   0.70711   0.00000
L2=   0.23570   0.23570   0.94281
L3=   0.66667   0.66667  -0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000
  0.00000000  0.04000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000 -0.04000000  0.00000000
  0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.18000001
TLS matrices from Libration in the M-base
  0.16722256  0.15722257 -0.07111052    0.04722228  0.03722228 -0.01111102    0.08277801  0.07277802 -0.04888849
  0.15722257  0.16722256 -0.07111052    0.03722228  0.04722228 -0.01111102    0.07277802  0.08277801 -0.04888849
 -0.07111052 -0.07111052  0.07555489   -0.01111102 -0.01111102  0.04555546   -0.04888849 -0.04888849 -0.01555602

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.16722256  0.15722257 -0.07111052    0.04722228  0.03722228 -0.01111102    0.08277801  0.07277802 -0.04888849
  0.15722257  0.16722256 -0.07111052    0.03722228  0.04722228 -0.01111102    0.07277802  0.08277801 -0.04888849
 -0.07111052 -0.07111052  0.07555489   -0.01111102 -0.01111102  0.04555546   -0.04888849 -0.04888849 -0.01555602

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=  -0.70711   0.70711   0.00000
Ly=   0.23570   0.23570   0.94281
Lz=   0.66667   0.66667  -0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.00999999  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000
  0.00000000  0.04000001  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000 -0.04000000  0.00000000
 -0.00000001  0.00000000  0.36000001    0.00000000  0.00000000  0.09000000    0.00000000  0.00000000  0.18000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.00000   0.00000   0.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    0.00000   0.00000   0.00000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   0.00000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  1.000000 -1.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.00000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000    0.01000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.04000000  0.00000000    0.00000000  0.04000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.36000001    0.00000000  0.00000000  0.36000001

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000   -0.00000001  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00000001  0.00000000  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000001 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000001

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000001 -0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000000  0.00000001
"""

getTLS_test194 = """\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
  -0.70711   0.70711   0.00000
   0.23570   0.23570   0.94281

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    1.00000   1.20000  -2.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   1.30000   0.50000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   1.50000
correlation shifts sx,sy,sz for the libration axes
  1.000000 -1.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.0000000 0.0000000 0.0000000
tx2,ty2,tz2= 0.0000000 0.0000000 0.0000000
vectors defining the principal Vibration axes
   1.00000   0.00000   0.00000
   0.00000   1.00000   0.00000

principal Libration axes (orthonormal L base)
L1=  -0.70711   0.70711   0.00000
L2=   0.23570   0.23570   0.94281
L3=   0.66667   0.66667  -0.33333

principal Vibration axes (orthonormal V base)
V1=   1.00000   0.00000   0.00000
V2=   0.00000   1.00000   0.00000
V3=   0.00000   0.00000   1.00000

TLS matrices from Libration in the L-base
  0.19004491  0.00000000 -0.00216858    0.01000000  0.00000000  0.00000000    0.01000000  0.02133329  0.01367080
  0.00000000  0.08551092  0.02916432    0.00000000  0.04000000  0.00000000   -0.00133340 -0.04000000  0.06505383
 -0.00216858  0.02916432  0.48448908    0.00000000  0.00000000  0.09000000    0.12727939  0.00000000  0.18000001
TLS matrices from Libration in the M-base
  0.32631192  0.13422243 -0.07313257    0.04722228  0.03722228 -0.01111102    0.02322223  0.13277806 -0.06499951
  0.13422246  0.32222280 -0.07211030    0.03722228  0.04722228 -0.01111102    0.03322224  0.16277806 -0.04299950
 -0.07313257 -0.07211030  0.11151022   -0.01111102 -0.01111102  0.04555546    0.02288919 -0.03888819 -0.03600030

V matrix from Vibration in the V-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
V matrix from Vibration in the M-base
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000
  0.00000000  0.00000000  0.00000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.32631192  0.13422243 -0.07313257    0.04722228  0.03722228 -0.01111102    0.02322223  0.13277806 -0.06499951
  0.13422246  0.32222280 -0.07211030    0.03722228  0.04722228 -0.01111102    0.03322224  0.16277806 -0.04299950
 -0.07313257 -0.07211030  0.11151022   -0.01111102 -0.01111102  0.04555546    0.02288919 -0.03888819 -0.03600030

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=  -0.70711   0.70711   0.00000
Ly=   0.23570   0.23570   0.94281
Lz=   0.66667   0.66667  -0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.19004489  0.00000000 -0.00216857    0.01000000  0.00000000  0.00000000    0.01000000  0.02133329  0.01367080
  0.00000000  0.08551093  0.02916432    0.00000000  0.04000000  0.00000000   -0.00133340 -0.04000000  0.06505383
 -0.00216860  0.02916433  0.48448908    0.00000000  0.00000000  0.09000000    0.12727939  0.00000000  0.18000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.14142  -1.36708   2.13333
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.62635   0.54212   0.03334
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   1.41422  -0.50000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    1.00000   1.20000  -2.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   1.30000   0.50000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  1.000000 -1.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.18004490  0.00000000 -0.00216858    0.01000000  0.00000000  0.00000000    0.19004491  0.00000000 -0.00216858
  0.00000000  0.04551093  0.02916432    0.00000000  0.04000000  0.00000000    0.00000000  0.08551092  0.02916432
 -0.00216858  0.02916432  0.12448908    0.00000000  0.00000000  0.36000001   -0.00216858  0.02916432  0.48448908

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.00000000  0.00000000  0.00000000   -0.00000001  0.00000000  0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000
  0.00000000  0.00000000  0.00000000   -0.00000002  0.00000001  0.00000000

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000001
  0.00000000  0.00000000  0.00000000    0.00000003 -0.00000001  0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   1.00000   0.00000   0.00000
Vy=   0.00000   1.00000   0.00000
Vz=   0.00000   0.00000   1.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.00000000  0.00000000  0.00000000    0.00000000 -0.00000001  0.00000001
  0.00000000  0.00000000  0.00000000    0.00000003 -0.00000001  0.00000001
  0.00000000  0.00000000  0.00000000    0.00000000  0.00000001  0.00000000
"""

getTLS_test195="""\n
** control information -NOT FOR COMPARISON- SKIP IT**

rms and rms2 Libration around i,j,k
dx ,dy ,dz = 0.1000000 0.2000000 0.3000000
dx2,dy2,dz2= 0.0100000 0.0400000 0.0900000
vectors defining the principal Libration axes
  -0.70711   0.70711   0.00000
   0.23570   0.23570   0.94281

rotation axes pass through the points in the M-system
 parallel to i : (wix,wiy,wiz)=    1.00000   1.20000  -2.00000
 parallel to j : (wjx,wjy,wjz)=   -1.00000   1.30000   0.50000
 parallel to k : (wkx,wky,wkz)=    0.00000   0.00000   1.50000
correlation shifts sx,sy,sz for the libration axes
  1.000000 -1.000000  2.000000

rms and rms2 Vibration along x,y,z
tx ,ty ,tz = 0.2000000 0.4000000 0.8000000
tx2,ty2,tz2= 0.0400000 0.1600000 0.6400000
vectors defining the principal Vibration axes
   0.00000   0.00000   1.00000
   1.00000   1.00000  -0.50000

principal Libration axes (orthonormal L base)
L1=  -0.70711   0.70711   0.00000
L2=   0.23570   0.23570   0.94281
L3=   0.66667   0.66667  -0.33333

principal Vibration axes (orthonormal V base)
V1=   0.00000   0.00000   1.00000
V2=   0.70711   0.70711   0.00000
V3=  -0.70711   0.70711   0.00000

TLS matrices from Libration in the L-base
  0.19004491  0.00000000 -0.00216858    0.01000000  0.00000000  0.00000000    0.01000000  0.02133329  0.01367080
  0.00000000  0.08551092  0.02916432    0.00000000  0.04000000  0.00000000   -0.00133340 -0.04000000  0.06505383
 -0.00216858  0.02916432  0.48448908    0.00000000  0.00000000  0.09000000    0.12727939  0.00000000  0.18000001
TLS matrices from Libration in the M-base
  0.32631192  0.13422243 -0.07313257    0.04722228  0.03722228 -0.01111102    0.02322223  0.13277806 -0.06499951
  0.13422246  0.32222280 -0.07211030    0.03722228  0.04722228 -0.01111102    0.03322224  0.16277806 -0.04299950
 -0.07313257 -0.07211030  0.11151022   -0.01111102 -0.01111102  0.04555546    0.02288919 -0.03888819 -0.03600030

V matrix from Vibration in the V-base
  0.04000000  0.00000000  0.00000000
  0.00000000  0.16000001  0.00000000
  0.00000000  0.00000000  0.64000005
V matrix from Vibration in the M-base
  0.40000001 -0.24000001  0.00000000
 -0.24000001  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000


************* INFORMATION FOR COMPARISON **********

***  T[M] L[M] S[M] *** total TLS matrices in the main base (initial information) ***
  0.72631192 -0.10577758 -0.07313257    0.04722228  0.03722228 -0.01111102    0.02322223  0.13277806 -0.06499951
 -0.10577755  0.72222281 -0.07211030    0.03722228  0.04722228 -0.01111102    0.03322224  0.16277806 -0.04299950
 -0.07313257 -0.07211030  0.15151022   -0.01111102 -0.01111102  0.04555546    0.02288919 -0.03888819 -0.03600030

***  Lx Ly Lz       *** principal Libration axes (orthonormal L base)
Lx=  -0.70711   0.70711   0.00000
Ly=   0.23570   0.23570   0.94281
Lz=   0.66667   0.66667  -0.33333

***  T[L] L[L] S[L] *** total TLS matrices in the L base ***
  0.83004487 -0.00000001 -0.00216855    0.01000000  0.00000000  0.00000000    0.01000000  0.02133329  0.01367080
  0.00000001  0.13884401  0.06687638    0.00000000  0.04000000  0.00000000   -0.00133340 -0.04000000  0.06505383
 -0.00216862  0.06687639  0.63115603    0.00000000  0.00000000  0.09000000    0.12727939  0.00000000  0.18000001

***  dx2 dy2 dz2     *** rms^2: Libration around lx,ly,lz
   0.0100000   0.0400000   0.0900000

***  dx  dy  dz      *** rms  : Libration around lx,ly,lz
   0.1000000   0.2000000   0.3000000

***  Wlx[L] Wly[L] Wlz[L] *** rotation axes pass through the points in the L-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    0.14142  -1.36708   2.13333
 Wly, axis parallel to ly : (wjx,wjy,wjz)=    1.62635   0.54212   0.03334
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   1.41422  -0.50000

***  Wlx[M] Wly[M] Wlz[M] *** rotation axes pass through the points in the M-base
 Wlx, axis parallel to lx : (wix,wiy,wiz)=    1.00000   1.20000  -2.00000
 Wly, axis parallel to ly : (wjx,wjy,wjz)=   -1.00000   1.30000   0.50000
 Wlz, axis parallel to lz : (wkx,wky,wkz)=    0.00000   0.00000   1.50000

***  sx sy sz              *** correlation shifts sx,sy,sz for the libration axes
  1.000000 -1.000000  2.000000

***  CW[L] CS[L] C[L]=CW[L]+CS[L] *** translation matrices from libration in the L base ***
  0.18004490  0.00000000 -0.00216858    0.01000000  0.00000000  0.00000000    0.19004491  0.00000000 -0.00216858
  0.00000000  0.04551093  0.02916432    0.00000000  0.04000000  0.00000000    0.00000000  0.08551092  0.02916432
 -0.00216858  0.02916432  0.12448908    0.00000000  0.00000000  0.36000001   -0.00216858  0.02916432  0.48448908

***  V[L]        V[L]=T[L]-C[L]     *** vibration matrix in the L-base
  0.63999999  0.00000000  0.00000001    0.63999999 -0.00000001  0.00000002
  0.00000000  0.05333309  0.03771205    0.00000001  0.05333309  0.03771206
  0.00000001  0.03771205  0.14666691   -0.00000004  0.03771207  0.14666694

***  V[M]        V[M]=RML*VM*RMLtr     *** vibration matrix in the M-base
  0.40000001 -0.24000001  0.00000000    0.40000001 -0.23999999  0.00000001
 -0.24000001  0.40000001  0.00000000   -0.23999994  0.40000001  0.00000000
  0.00000000  0.00000000  0.04000000   -0.00000001  0.00000001  0.04000000

***  Vx Vy Vz       *** principal Vibration axes (orthonormal V base)
Vx=   0.00000   0.00000   1.00000
Vy=   0.70711   0.70711   0.00000
Vz=  -0.70711   0.70711   0.00000

***  V[V]        V[V]=RMVtr*VM*RMV  *** vibration matrix in the V-base
  0.04000000  0.00000000  0.00000000    0.04000000  0.00000000  0.00000002
  0.00000000  0.16000001  0.00000000    0.00000001  0.16000006 -0.00000004
  0.00000000  0.00000000  0.64000005   -0.00000001  0.00000002  0.63999993
"""

def generate_input_pdb(T,L,S, pdb_lines, prefix):
  pi = iotbx.pdb.input(source_info=None, lines=pdb_lines)
  ph = pi.construct_hierarchy()
  xrs = pi.xray_structure_simple()
  sel = [flex.bool(xrs.scatterers().size(), True).iselection()]
  tlsos = mmtbx.tls.tools.generate_tlsos(
    selections     = sel,
    xray_structure = xrs,
    T=[T.as_sym_mat3()],
    L=[L.as_sym_mat3()],
    S=[S])
  u_cart_from_tls = mmtbx.tls.tools.u_cart_from_tls(
    sites_cart = xrs.sites_cart(),
    selections = sel,
    tlsos      = tlsos)
  xrs.convert_to_anisotropic()
  xrs.set_u_cart(u_cart=u_cart_from_tls)
  ph.adopt_xray_structure(xrs)
  of = open("%s_start.pdb"%prefix, "w")
  mmtbx.tls.tools.remark_3_tls(tlsos=tlsos, selection_strings=["all"], out=of)
  of.close()
  ph.write_pdb_file(
    file_name        = "%s_start.pdb"%prefix,
    open_append      = True,
    crystal_symmetry = xrs.crystal_symmetry())

def generate_ensemble(pdb_lines, r, prefix):
  # read inputs
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_lines)
  xrs = pdb_inp.xray_structure_simple()
  xrs.convert_to_isotropic()
  xrs = xrs.set_b_iso(value=0)
  sites_cart = xrs.sites_cart()
  h = pdb_inp.construct_hierarchy()
  #states = mmtbx.utils.states(xray_structure=xrs, pdb_hierarchy=h)
  #
  #remark_3_records = pdb_inp.extract_remark_iii_records(3)
  #tls_extract = mmtbx.tls.tools.tls_from_pdb_inp(
  #  remark_3_records = remark_3_records,
  #  pdb_hierarchy = h)
  #assert len(tls_extract.tls_params)==1
  #tlso = tls_extract.tls_params[0]
  #print "Input T:", tlso.t
  #print "      L:", tlso.l
  #print "      S:", tlso.s
  #T = matrix.sym(sym_mat3=tlso.t)
  #L = matrix.sym(sym_mat3=tlso.l)#*sc #???
  #S = matrix.sqr(tlso.s)
  # run through steps
  #r = tls_as_xyz.decompose_tls(T=T, L=L, S=S)
  # iteratable
  tls_as_xyz.ensemble_generator(
    decompose_tls_object = r,
    pdb_hierarchy        = h,
    xray_structure       = xrs,
    n_models             = 499, # that's how many my version of Pymol can show
    log                  = sys.stdout).write_pdb_file(
      file_name="%s_ensemble.pdb"%prefix)

def run_core(lines, prefix, comp, generate):
  log = open("%s.log"%prefix, "w")
  print_step("Test %s:"%prefix, log)
  e = extract(s=lines)
  r = tls_as_xyz.decompose_tls(T=e.T_M, L=e.L_M, S=e.S_M, log=log)
  log.close()
  if(comp): compare(e, r)
  if(generate):
    generate_input_pdb(T=e.T_M,L=e.L_M,S=e.S_M,pdb_lines=pdb_str,prefix=prefix)
    generate_ensemble(pdb_lines=pdb_str, r=r, prefix=prefix)

def exercises_0():
  for lines, prefix in [(getTLS3_test001, "getTLS3_test001"),
                        (getTLS3_test004, "getTLS3_test004"),
                        (getTLS3_test011, "getTLS3_test011"),
                        (getTLS3_test014, "getTLS3_test014"),
                        (getTLS3_test021, "getTLS3_test021"),
                        (getTLS3_test024, "getTLS3_test024"),
                        (getTLS3_test031, "getTLS3_test031"),
                        (getTLS3_test034, "getTLS3_test034"),
                        (getTLS3_test041, "getTLS3_test041"),
                        (getTLS3_test044, "getTLS3_test044")]:
    run_core(lines=lines, prefix=prefix, comp=True, generate=False)

def exercises_1():
  for lines, prefix in [(getTLS_test111, "getTLS_test111"),
                        (getTLS_test112, "getTLS_test112"),
                        (getTLS_test113, "getTLS_test113"),
                        (getTLS_test114, "getTLS_test114"),
                        (getTLS_test121, "getTLS_test121"),
                        (getTLS_test122, "getTLS_test122"),
                        (getTLS_test123, "getTLS_test123"),
                        (getTLS_test124, "getTLS_test124"),
                        (getTLS_test131, "getTLS_test131"),
                        (getTLS_test132, "getTLS_test132"),
                        (getTLS_test133, "getTLS_test133"),
                        (getTLS_test134, "getTLS_test134"),

                        (getTLS_test141, "getTLS_test141"),
                        (getTLS_test142, "getTLS_test142"),
                        (getTLS_test143, "getTLS_test143"),
                        (getTLS_test144, "getTLS_test144")
                       ]:
    run_core(lines=lines, prefix=prefix, comp=True, generate=True)

def exercises_2():
  # ensembles from 16x series must be exactly the same as from 17x
  for lines, prefix in [(getTLS_test161, "getTLS_test161"),
                        (getTLS_test162, "getTLS_test162"),
                        (getTLS_test163, "getTLS_test163"),
                        (getTLS_test164, "getTLS_test164")
                       ]:
    run_core(lines=lines, prefix=prefix, comp=False, generate=True)
  for lines, prefix in [(getTLS_test171,    "getTLS_test171"),
                        (getTLS_test172,    "getTLS_test172"),
                        (getTLS_test173,    "getTLS_test173"),
                        (getTLS_test174,    "getTLS_test174"),
                        (getTLS_test174new, "getTLS_test174new"),
                        (getTLS_test175,    "getTLS_test175"),
                        (getTLS_test176,    "getTLS_test176")]:
    run_core(lines=lines, prefix=prefix, comp=False, generate=True)
    # XXX Also check ensembles 16x == 17x

def exercises_3():
  # Vot seriya 18x; napravleniya vekrotov translyacii i vrascheniya odinakovy
  # dlya vseh treh sluchaev, no v 181 vraschenie >> translyacii, v 182
  # vreaschenie << translyacii, v 183 oni sravnimy (kam mne kazhetsya
  # teoreticheski).
  # Da, kazhetsya , my esche ne probovali vraschenie srazu vokrug treh osej ?
  # Dlya etogo seriya 19x. No tut glazom uzhe slozhno proveryatx. Razve chto
  # sravnitx tvoi i moj chisla ?
  for lines, prefix in [(getTLS_test181,    "getTLS_test181"),
                        (getTLS_test181new, "getTLS_test181new"),
                        (getTLS_test182,    "getTLS_test182"),
                        (getTLS_test183,    "getTLS_test183"),

                        (getTLS_test191,    "getTLS_test191"),
                        #(getTLS_test192,    "getTLS_test192"),
                        #(getTLS_test193,    "getTLS_test193"),
                        #(getTLS_test194,    "getTLS_test194"),
                        #(getTLS_test195,    "getTLS_test195")
                       ]:
    run_core(lines=lines, prefix=prefix, comp=True, generate=True)

if (__name__ == "__main__"):
  if(1): exercises_0()
  if(0): exercises_1() # currently fails
  if(1): exercises_2()
  if(1): exercises_3()
