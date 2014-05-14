from __future__ import division
import sys

from mmtbx.conformation_dependent_library.cdl_database import cdl_database

step = 10

def round_to_ten(d):
  t = int(d//10)*10
  if t==180: return -180
  return t

def get_grid_values(phi, psi, column=2):
  print 'phi',phi
  print 'psi',psi
  key0 = (round_to_ten(phi), round_to_ten(psi))
  grid = []
  indices = []
  for j in range(-1,3):
    grid.append([])
    indices.append([])
    for i in range(-1,3):
      key = (key0[0]+i*step, key0[1]+j*step)
      grid[-1].append(cdl_database["Gly_nonxpro"][key][column])
      indices[-1].append(key)
  if 1:
    for i, d in zip(indices, grid):
      print i,d
  return grid

def get_index(phi, psi):
  key0 = (round_to_ten(phi), round_to_ten(psi))
  index = ((phi-key0[0])/step+1, (psi-key0[1])/10+1)
  return index


def run():
  print 'running'
  for i in range(91,120):
    for j in range(1,20):
      grid = get_grid_values(float(i),float(j))
      print grid
      index = get_index(float(i), float(j))
      print index
      assert 0

# XXX duplicates scitbx/math/interpolation.h
def interpolate_at_point (p0, p1, p2, p3, t) :
  """
  http://en.wikipedia.org/wiki/Cubic_Hermite_spline
  http://www.mvps.org/directx/articles/catmull
  """
  t2 = t * t
  t3 = t2 * t
  result = 0.5 * ((2.0*p1) + (-p0 + p2) * t + \
       (2.0*p0 - 5.0*p1 + 4.0*p2 - p3) * t2 + \
       (-p0 + 3.0*p1 - 3.0*p2 + p3) * t3)
  return result

def interpolate_2d (values, coords) :
  """
  Given a 2-dimensional (4x4) grid of values and fractional coordinates,
  compute the value at these coordinates.  Assumes that the grid covers the
  points from -1 to 2 in both dimensions, and the fracional coordinates are
  between 0 and 1.
  """
  x, y = coords
  assert (0 <= x <= 1) and (0 <= y <= 1)
  p = []
  for i in range(0, 4) :
    p0, p1, p2, p3 = values[i]
    p.append(interpolate_at_point(p0, p1, p2, p3, x))
  p0, p1, p2, p3 = p
  result = interpolate_at_point(p0, p1, p2, p3, y)
  return result

def exercise () :
  values = [
    [0, 1, 2, 1],
    [0, 2, 2, 1.5],
    [1, 1, 3, 1],
    [1, 2, 3, 0] ]
  x, y = (0.5, 0.5)
  r = interpolate_2d(values, (x, y))
  print r

if __name__=="__main__":
  #exercise()
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
