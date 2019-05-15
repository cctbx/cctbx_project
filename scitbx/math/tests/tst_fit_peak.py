from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex
from scitbx.math import fit_peak
from libtbx.test_utils import approx_equal
from six.moves import range

# =============================================================================
def test_pick_map_neighbors():

  # generate test map
  d = 3
  gridding = flex.grid(d,d,d)
  test_map = flex.double(gridding)
  count = 1
  for i in range(d):
    for j in range(d):
      for k in range(d):
        test_map[(i,j,k)] = count
        count = count + 1

  # order neighbors
  fitted = fit_peak.pick_map_neighbors(site=(1,1,1),map_in=test_map)
  height_list,xyz_list = fitted.get_26_nearest_neighbors()

  # check output
  height_list_correct = [14.0,23.0,5.0,17.0,11.0,15.0,13.0,18.0,4.0,20.0,12.0,
                         22.0,2.0,16.0,6.0,26.0,10.0,24.0,8.0,27.0,9.0,21.0,
                         25.0,3.0,7.0,19.0,1.0]
  xyz_list_correct = [(1.0/3.0,1.0/3.0,1.0/3.0),(2.0/3.0,1.0/3.0,1.0/3.0),
                      (0.0,1.0/3.0,1.0/3.0),(1.0/3.0,2.0/3.0,1.0/3.0),
                      (1.0/3.0,0.0,1.0/3.0),(1.0/3.0,1.0/3.0,2.0/3.0),
                      (1.0/3.0,1.0/3.0,0.0),(1.0/3.0,2.0/3.0,2.0/3.0),
                      (0.0,1.0/3.0,0.0),(2.0/3.0,0.0,1.0/3.0),
                      (1.0/3.0,0.0,2.0/3.0),(2.0/3.0,1.0/3.0,0.0),
                      (0.0,0.0,1.0/3.0),(1.0/3.0,2.0/3.0,0.0),
                      (0.0,1.0/3.0,2.0/3.0),(2.0/3.0,2.0/3.0,1.0/3.0),
                      (1.0/3.0,0.0,0.0),(2.0/3.0,1.0/3.0,2.0/3.0),
                      (0.0,2.0/3.0,1.0/3.0),(2.0/3.0,2.0/3.0,2.0/3.0),
                      (0.0,2.0/3.0,2.0/3.0),(2.0/3.0,0.0,2.0/3.0),
                      (2.0/3.0,2.0/3.0,0.0),(0.0,0.0,2.0/3.0),
                      (0.0,2.0/3.0,0.0),(2.0/3.0,0.0,0.0),(0.0,0.0,0.0)]
  for i in range(27):
    assert(height_list[i] == height_list_correct[i])
    assert(xyz_list[i] == xyz_list_correct[i])

# =============================================================================
def test_fit_3d_parabola():

  # generate test map
  d = 10
  gridding = flex.grid(d,d,d)
  test_map = flex.double(gridding)
  parameters = (1.5,1.5,1.5,-2.0,-2.0,-2.0,24.0)
  p = fit_peak.parabola(parameters=parameters)
  for i in range(d):
    for j in range(d):
      for k in range(d):
        test_map[(i,j,k)] = p.get_height(r=(float(i)/d,float(j)/d,float(k)/d))

  # check fitting
  neighbors = fit_peak.pick_map_neighbors(site=(6,6,6),map_in=test_map)
  height_list,xyz_list = neighbors.get_26_nearest_neighbors()
  fp = fit_peak.fit_peak(height_list=height_list,xyz_list=xyz_list,
                         shape="parabola")
  assert approx_equal(fp.x, parameters)
  assert approx_equal(fp.vertex, p.vertex)

# =============================================================================
def test_fit_3d_quadratic():

  # generate test map
  d = 10
  gridding = flex.grid(d,d,d)
  test_map = flex.double(gridding)
  parameters = (1.5,1.5,1.5,-2.0,-2.0,-2.0,1.25,1.25,1.25,24.0)
  p = fit_peak.quadratic(parameters=parameters)
  for i in range(d):
    for j in range(d):
      for k in range(d):
        test_map[(i,j,k)] = p.get_height(r=(float(i)/d,float(j)/d,float(k)/d))

  # check fitting
  neighbors = fit_peak.pick_map_neighbors(site=(4,4,4),map_in=test_map)
  height_list,xyz_list = neighbors.get_26_nearest_neighbors()
  fp = fit_peak.fit_peak(height_list=height_list,xyz_list=xyz_list,
                         shape="quadratic")
  assert approx_equal(fp.x, parameters)
  assert approx_equal(fp.vertex, p.vertex)

# =============================================================================
def test_fit_3d_gaussian():

  # generate test map
  d = 10
  gridding = flex.grid(d,d,d)
  test_map = flex.double(gridding)
  parameters = (1.5,1.5,1.5,-2.0,-2.0,-2.0,1.25,1.25,1.25,24.0)
  p = fit_peak.gaussian(parameters=parameters)
  for i in range(d):
    for j in range(d):
      for k in range(d):
        test_map[(i,j,k)] = p.get_height(r=(float(i)/d,float(j)/d,float(k)/d))

  # check fitting
  neighbors = fit_peak.pick_map_neighbors(site=(4,4,4),map_in=test_map)
  height_list,xyz_list = neighbors.get_26_nearest_neighbors()
  fp = fit_peak.fit_peak(height_list=height_list,xyz_list=xyz_list,
                         shape="gaussian")
  assert approx_equal(fp.x, parameters)
  assert approx_equal(fp.vertex, p.vertex)

# =============================================================================
if (__name__ == "__main__"):
  test_pick_map_neighbors()
  test_fit_3d_parabola()
  test_fit_3d_quadratic()
  test_fit_3d_gaussian()
  print("OK")
