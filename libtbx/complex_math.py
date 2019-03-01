from __future__ import absolute_import, division, print_function
from math import pi, atan2, cos, sin

def norm(c):
  return c.real**2 + c.imag**2

def abs_arg(c, deg=False):
  "conversion of complex number: real, imag -> absolute value, polar angle"
  a = abs(c)
  if (a == 0): return (0, 0)
  t = atan2(c.imag, c.real)
  if (deg): t *= 180/pi
  return (a, t)

def arg(c, deg=False):
  "conversion of complex number: real, imag -> polar angle"
  return abs_arg(c, deg)[1]

def polar(a_t, deg=False):
  "conversion of complex number: polar representation -> real, imag"
  a, t = a_t
  if (deg): t *= pi/180
  return complex(a * cos(t), a * sin(t))
