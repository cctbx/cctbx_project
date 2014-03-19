from __future__ import division
from scitbx.array_family import flex
import iotbx.pdb
from cctbx import maptbx
from cctbx import miller

def getvs(cmap, threshold):
  co = maptbx.connectivity(map_data=cmap, threshold=threshold)
  map_result = co.result()
  v=[0,0,0]
  for i in range(3):
    v[i] = (map_result==i).count(True)
  return v


def exercise1():
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       2.000   2.000   2.000  1.00 20.00           C
END
"""
  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  cg = maptbx.crystal_gridding(unit_cell=xrs.unit_cell(),
      pre_determined_n_real=(100,100,100),
      space_group_info=xrs.space_group_info())
  fc = xrs.structure_factors(d_min = 1., algorithm = "direct").f_calc()
  fft_map = miller.fft_map(crystal_gridding=cg, fourier_coefficients=fc)
  map_data = fft_map.real_map_unpadded()
  v = getvs(map_data, 100)
  assert v[2] == 0
  assert v[1] < 15000
  assert v[0]+v[1]+v[2] == 1000000

def exercise3():
  pdb_str="""
CRYST1   10.000  10.000   10.000  90.00  90.00  90.00 P 1
HETATM    1  C    C      1       2.000   2.000   2.000  1.00  2.00           C
HETATM    1  C    C      1       3.500   2.000   2.000  1.00  2.00           C
END
"""

  pdb_inp = iotbx.pdb.input(source_info=None, lines=pdb_str)
  xrs = pdb_inp.xray_structure_simple()
  cg = maptbx.crystal_gridding(unit_cell=xrs.unit_cell(),
      pre_determined_n_real=(100,100,100),
      space_group_info=xrs.space_group_info())
  fc = xrs.structure_factors(d_min = 1., algorithm = "direct").f_calc()
  fft_map = miller.fft_map(crystal_gridding=cg, fourier_coefficients=fc)
  fft_map.apply_sigma_scaling()
  map_data = fft_map.real_map_unpadded()
  #all filled
  v = getvs(map_data, -100)
  assert v == [0, 1000000, 0]
  # can see one blob
  v = getvs(map_data, 5)
  assert v[0]+v[1]+v[2] == 1000000
  assert v[2] == 0
  # can see separate, approx equal volume bloobs
  v = getvs(map_data, 10)
  assert v[0]+v[1]+v[2] == 1000000
  assert abs(v[1] - v[2]) < 5
  # nothing to see
  v = getvs(map_data, 1000)
  assert v == [1000000, 0, 0]

def exercise4():
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(1)
  for i in range(10,20):
    for j in range(10,20):
      for k in range(10,20):
        cmap[i,j,k] = 10
  v = getvs(cmap, 5)
  assert v == [999000, 1000, 0]
  #print "all filled"
  v = getvs(cmap, -5)
  assert v == [0,1000000,0]
  #print "none filled"
  v = getvs(cmap, 20)
  assert v == [1000000,0,0]

def exercise5():
  #print "corner blob"
  cmap = flex.double(flex.grid(100,100,100))
  cmap.fill(0)
  for i in range(100):
    for j in range(100):
      for k in range(100):
        if (i<10 or i>=90) and (j<10 or j>=90) and (k<10 or k>=90):
          cmap[i,j,k] = 10
          #print i,j,k
  v = getvs(cmap, 5)
  assert v == [992000, 8000, 0]

  #print "2 blobs"
  cmap.fill(0)
  for i in range(100):
    for j in range(100):
      for k in range(100):
        if (5<i<10) and (5<j<10) and (5<k<10):
          cmap[i,j,k] = 10
        if (15<i<20) and (15<j<20) and (15<k<20):
          cmap[i,j,k] = 20
  v = getvs(cmap, 5)
  assert v == [999872,64,64]
  v = getvs(cmap, 15)
  assert v == [999936, 64,0]

  #print "endless blob"
  cmap.fill(0)
  for j in range(100):
    for k in range(100):
      cmap[5,j,k] = 10
  v = getvs(cmap, 5)
  assert v == [990000, 10000, 0]

def exercise6():
  cmap = flex.double(flex.grid(100,100,100))
  #print "corner touch"
  cmap.fill(0)
  cmap[1,1,1] = cmap[2,2,2] = 10
  v = getvs(cmap, 5)
  assert v == [999998, 1, 1]
  #print "edges touch"
  cmap.fill(0)
  cmap[1,1,1] = cmap[2,2,1] = 10
  v = getvs(cmap, 5)
  assert v == [999998, 1, 1]
  #print "face touch"
  cmap.fill(0)
  cmap[1,1,1] = cmap[2,1,1] = 10
  v = getvs(cmap, 5)
  assert v == [999998, 2, 0]


if __name__ == "__main__" :
  exercise1()
  exercise3()
  exercise4()
  exercise5()
  exercise6()
