from math import pi, atan2, cos, sin

def abs_arg(c, deg=False):
  a = abs(c)
  if (a == 0): return (0, 0)
  t = atan2(c.imag, c.real)
  if (deg): t *= 180/pi
  return (a, t)

def polar(a_t, deg=False):
  a, t = a_t
  if (deg): t *= pi/180
  return complex(a * cos(t), a * sin(t))
