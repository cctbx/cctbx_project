from scitbx.math import euler_angles_xyz_matrix, euler_angles_xyz_angles
from scitbx.math import euler_angles_yzx_matrix, euler_angles_yzx_angles
from scitbx.math import euler_angles_zyz_matrix, euler_angles_zyz_angles

"""
Mathematica code:
rx={{1,0,0},{0,cx,-sx},{0,sx,cx}}
ry={{cy,0,sy},{0,1,0},{-sy,0,cy}}
rz={{cz,-sz,0},{sz,cz,0},{0,0,1}}
rz1={{cz1,-sz1,0},{sz1,cz1,0},{0,0,1}}
rz3={{cz3,-sz3,0},{sz3,cz3,0},{0,0,1}}
"""

def xyz_matrix(ax, ay, az):
  """
  2008: C++ implementation 
  """
  return euler_angles_xyz_matrix( ax, ay, az )

def xyz_angles(m, eps=1.e-12):
  """
  2008: C++ implementation 
  """
  return euler_angles_xyz_angles( m, eps )

def yzx_matrix(ay, az, ax):
  """
  2008: C++ implementation 
  """
  return euler_angles_yzx_matrix( ay, az, ax )

def yzx_angles(m, eps=1.e-12):
  """
  2008: C++ implementation 
  """
  return euler_angles_yzx_angles( m, eps )

def zyz_matrix(az1, ay2, az3):
  """
  2008: C++ implementation 
  """
  return euler_angles_zyz_matrix( az1, ay2, az3 )

def zyz_angles(m, eps=1.e-12):
  """
  2008: C++ implementation 
  """
  return euler_angles_zyz_angles( m, eps )
