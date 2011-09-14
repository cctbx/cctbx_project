from cctbx.sgtbx.direct_space_asu import direct_space_asu
from cctbx.sgtbx.direct_space_asu.short_cuts import *

def asu_01(): # p_1 (s.g. 1)
  return (direct_space_asu('P 1')
    & x0
    & +x1
    & y0
    & +y1
    & z0
    & +z1
)

def asu_02(): # p_2 (s.g. 3)
  return (direct_space_asu('P 2')
    & x0(y2)
    & x2(y2)
    & y0
    & +y1
    & z0
    & +z1
)

def asu_03(): # p_m (s.g. 6)
  return (direct_space_asu('P -2x')
    & x0
    & x2
    & y0
    & +y1
    & z0
    & +z1
)

def asu_04(): # p_g (s.g. 7)
  return (direct_space_asu('P -2xb')
    & x0(+y2)
    & x2(+y2)
    & y0
    & +y1
    & z0
    & +z1
)

def asu_05(): # c_m (s.g. 8)
  return (direct_space_asu('C -2x')
    & x0
    & x2
    & y0
    & +y2
    & z0
    & +z1
)

def asu_06(): # p_2_m_m (s.g. 25)
  return (direct_space_asu('P 2 -2')
    & x0
    & x2
    & y0
    & y2
    & z0
    & +z1
)

def asu_07(): # p_2_m_g (s.g. 28)
  return (direct_space_asu('P 2 -2a')
    & x0(y2)
    & x4
    & y0
    & +y1
    & z0
    & +z1
)

def asu_08(): # p_2_g_g (s.g. 32)
  return (direct_space_asu('P 2 -2ab')
    & x0
    & x2(-y0)
    & y0
    & +y2
    & z0
    & +z1
)

def asu_09(): # c_2_m_m (s.g. 35)
  return (direct_space_asu('C 2 -2')
    & x0
    & x4(y4)
    & y0
    & y2
    & z0
    & +z1
)

def asu_10(): # p_4 (s.g. 75)
  return (direct_space_asu('P 4')
    & x0(-y0)
    & x2
    & y0
    & y2(-x2)
    & z0
    & +z1
)

def asu_11(): # p_4_m_m (s.g. 99)
  return (direct_space_asu('P 4 -2')
    & x0
    & y2
    & -p0
    & z0
    & +z1
)

def asu_12(): # p_4_g_m (s.g. 100)
  return (direct_space_asu('P 4 -2ab')
    & x0(-y0)
    & y0
    & m2
    & z0
    & +z1
)

def asu_13(): # p_3 (s.g. 143)
  return (direct_space_asu('P 3')
    & x0(-y0)
    & y0
    & k1
    & m1(-h1 | -k1)
    & h1
    & z0
    & +z1
)

def asu_14(): # p_3_m_1 (s.g. 156)
  return (direct_space_asu('P 3 -2"')
    & h0
    & m1
    & k0
    & z0
    & +z1
)

def asu_15(): # p_3_1_m (s.g. 157)
  return (direct_space_asu('P 3 -2')
    & y0
    & k1
    & m1(y3)
    & p0
    & z0
    & +z1
)

def asu_16(): # p_6 (s.g. 168)
  return (direct_space_asu('P 6')
    & y0
    & k1
    & m1(y3)
    & p0(-y0)
    & z0
    & +z1
)

def asu_17(): # p_6_m_m (s.g. 183)
  return (direct_space_asu('P 6 -2')
    & y0
    & k1
    & -h0
    & z0
    & +z1
)

def get_asu(point_group_number):
  return eval("asu_%02d" % point_group_number)()

if (__name__ == "__main__"):
  for i in xrange(1,17+1):
    get_asu(i).show_summary()
