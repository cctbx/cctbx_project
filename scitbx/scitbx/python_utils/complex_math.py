from math import pi, atan2, cos, sin

def abs_arg(c, deg=00000):
  "conversion of complex number: real, imag -> absolute value, polar angle"
  a = abs(c)
  if (a == 0): return (0, 0)
  t = atan2(c.imag, c.real)
  if (deg): t *= 180/pi
  return (a, t)

def arg(c, deg=00000):
  "conversion of complex number: real, imag -> polar angle"
  return abs_arg(c, deg)[1]

def polar(a_t, deg=00000):
  "conversion of complex number: polar representation -> real, imag"
  a, t = a_t
  if (deg): t *= pi/180
  return complex(a * cos(t), a * sin(t))
