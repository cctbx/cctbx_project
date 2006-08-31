import math

"""
Mathematica code:
rx={{1,0,0},{0,cx,-sx},{0,sx,cx}}
ry={{cy,0,sy},{0,1,0},{-sy,0,cy}}
rz={{cz,-sz,0},{sz,cz,0},{0,0,1}}
rz1={{cz1,-sz1,0},{sz1,cz1,0},{0,0,1}}
rz3={{cz3,-sz3,0},{sz3,cz3,0},{0,0,1}}
"""

def xyz_matrix(ax, ay, az):
  ax *= math.pi/180
  ay *= math.pi/180
  az *= math.pi/180
  sx = math.sin(ax)
  sy = math.sin(ay)
  sz = math.sin(az)
  cx = math.cos(ax)
  cy = math.cos(ay)
  cz = math.cos(az)
  # rx.ry.rz
  return (
     cy*cz,             -cy*sz,              sy,
     cx*sz + sx*sy*cz,   cx*cz - sx*sy*sz,  -sx*cy,
     sx*sz - cx*sy*cz,   sx*cz + cx*sy*sz,   cx*cy)

def xyz_angles(m, eps=1.e-12):
  m00, m01, m02, \
  m10, m11, m12, \
  m20, m21, m22 = m
  if   (m02 >  1-eps):
    ax = math.atan2(m21,m11)
    ay = math.pi/2
    az = 0.
  elif (m02 < -1+eps):
    ax = math.atan2(m21,m11)
    ay = -math.pi/2
    az = 0.
  else:
    ax = math.atan2(-m12,m22)
    ay = math.asin(m02)
    az = math.atan2(-m01,m00)
  return [v*180/math.pi for v in ax, ay, az]

def yzx_matrix(ay, az, ax):
  ay *= math.pi/180
  az *= math.pi/180
  ax *= math.pi/180
  cy = math.cos(ay)
  cz = math.cos(az)
  cx = math.cos(ax)
  sy = math.sin(ay)
  sz = math.sin(az)
  sx = math.sin(ax)
  # ry.rz.rx
  return (
     cy*cz,   sx*sy - cx*cy*sz,   cx*sy + cy*sx*sz,
     sz,      cx*cz,             -cz*sx,
    -cz*sy,   cy*sx + cx*sy*sz,   cx*cy - sx*sy*sz)

def yzx_angles(m, eps=1.e-12):
  m00, m01, m02, \
  m10, m11, m12, \
  m20, m21, m22 = m
  if   (m10 >  1-eps):
    ay = math.atan2(m02,m22)
    az = math.pi/2
    ax = 0.
  elif (m10 < -1+eps):
    ay = math.atan2(m02,m22)
    az = -math.pi/2
    ax = 0.
  else:
    ay = math.atan2(-m20,m00)
    az = math.asin(m10)
    ax = math.atan2(-m12,m11)
  return [v*180/math.pi for v in ay, az, ax]

def zyz_matrix(az1, ay, az3):
  az1 *= math.pi/180
  ay  *= math.pi/180
  az3 *= math.pi/180
  sz1 = math.sin(az1)
  sy  = math.sin(ay)
  sz3 = math.sin(az3)
  cz1 = math.cos(az1)
  cy  = math.cos(ay)
  cz3 = math.cos(az3)
  # rz1.ry.rz3
  return (
     cy*cz1*cz3 - sz1*sz3,  -cz3*sz1 - cy*cz1*sz3,   cz1*sy,
     cy*cz3*sz1 + cz1*sz3,   cz1*cz3 - cy*sz1*sz3,   sy*sz1,
    -cz3*sy,                 sy*sz3,                 cy)

def zyz_angles(m, eps=1.e-12):
  m00, m01, m02, \
  m10, m11, m12, \
  m20, m21, m22 = m
  if   (m22 >  1-eps):
    az1 = math.atan2(-m01,m11)
    ay  = 0.
    az3 = 0.
  elif (m22 < -1+eps):
    az1 = math.atan2(-m01,m11)
    ay  = math.pi
    az3 = 0.
  else:
    az1 = math.atan2(m12,m02)
    ay  = math.acos(m22)
    az3 = math.atan2(m21,-m20)
  return [v*180/math.pi for v in az1, ay, az3]
