from __future__ import absolute_import, division, print_function
from mmtbx.conformation_dependent_library.cdl_database import cdl_database
from six.moves import range

step = 10

def round_to_ten(d):
  t = int(d//10)*10
  if t==180: return -180
  return t

def get_db_result(db, key, column):
  if key in db:
    return db[key][column]
  key = list(key)
  for i in range(2):
    if key[i]<-180: key[i]+=360
    if key[i]>170: key[i]-=360
  return db[tuple(key)][column]

def get_grid_values(residue_type, phi, psi, column=2):
  key0 = (round_to_ten(phi), round_to_ten(psi))
  grid = []
  indices = []
  for j in range(-1,3):
    grid.append([])
    indices.append([])
    for i in range(-1,3):
      key = (key0[0]+i*step, key0[1]+j*step)
      grid[-1].append(get_db_result(cdl_database[residue_type], key, column))
      indices[-1].append(key)
  return grid

def get_index(phi, psi):
  key0 = (round_to_ten(phi), round_to_ten(psi))
  if phi>=180: phi-=360
  if psi>=180: psi-=360
  index = ((phi-key0[0])/step,
           (psi-key0[1])/step,
          )
  return index

def print_grid(grid, phi, psi):
  outl = "-"*30
  outl += "\n"
  for i in range(4):
    for j in range(4):
      outl += " %6.2f" % grid[i][j]
    outl += "\n"
  outl += "-"*30
  outl += "\n"
  print(outl)

# XXX duplicates scitbx/math/interpolation.h
def interpolate_at_point(p0, p1, p2, p3, t):
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

def interpolate_2d(values, coords):
  """
  Given a 2-dimensional (4x4) grid of values and fractional coordinates,
  compute the value at these coordinates.  Assumes that the grid covers the
  points from -1 to 2 in both dimensions, and the fracional coordinates are
  between 0 and 1.
  """
  x, y = coords
  assert (0 <= x <= 1) and (0 <= y <= 1)
  p = []
  for i in range(0, 4):
    p0, p1, p2, p3 = values[i]
    p.append(interpolate_at_point(p0, p1, p2, p3, x))
  p0, p1, p2, p3 = p
  result = interpolate_at_point(p0, p1, p2, p3, y)
  return result
