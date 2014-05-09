from __future__ import division
import os, sys

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
  for j in range(-2,2):
    for i in range(-1,3):
      key = (key0[0]+i*step, key0[1]+j*step)
      grid.append(cdl_database["Gly_nonxpro"][key][column])
      indices.append(key)
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

if __name__=="__main__":
  args = sys.argv[1:]
  del sys.argv[1:]
  run(*tuple(args))
