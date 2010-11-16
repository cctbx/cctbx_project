import iotbx.xplor.map
from cctbx import maptbx
from cctbx import sgtbx
from cctbx import uctbx
from cctbx.array_family import flex
from cctbx.development import random_structure
from libtbx.math_utils import iround
from libtbx.test_utils import approx_equal, eps_eq, show_diff
import libtbx.load_env
from cStringIO import StringIO
import os

def exercise_map_gridding():
  try:
    g = iotbx.xplor.map.gridding(n=(0,20,30), first=(-3,-4,-5), last=(5,4,3))
  except RuntimeError, e:
    assert str(e) == "Illegal xplor map gridding for dimension X: " \
                   + "gridding=0, first=-3, last=5"
  g = iotbx.xplor.map.gridding(n=(10,20,30), first=(-3,-4,-5), last=(5,4,3))
  fg = g.as_flex_grid()
  assert fg.origin() == g.first
  assert fg.last(False) == g.last
  assert fg.last(True) == (6,5,4)
  assert fg.all() == (9,9,9)
  assert not fg.is_padded()

def get_test_file_name():
  """CNS commands used to generate test map:
xray
  MAPResolution 1
  declare name=map domain=real end
  do (map=ft(fcalc)) (all)
  write map
    from=map output=cns.map extend=box
    xmin=1 xmax=3
    ymin=-1 ymax=0
    zmin=-2 zmax=-1
  end
end
  """
  return libtbx.env.find_in_repositories(
    relative_path="phenix_regression/misc/cns.map", test=os.path.isfile)

def read_xplor(file_name):
  a = iotbx.xplor.map.reader(file_name=file_name)
  assert a.title_lines == [
    ' REMARKS FILENAME="cns.map"',
    ' REMARKS DATE:15-May-2004  02:15:56       created by user: rwgk',
    ' REMARKS VERSION:1.1']
  assert a.gridding.n == (24,24,40)
  assert a.gridding.first == (1,-4,-6)
  assert a.gridding.last == (10,0,-3)
  assert a.unit_cell.is_similar_to(
    uctbx.unit_cell((7.41407939496,7.41407939496,12.6039349714,90,90,120)))
  assert approx_equal(a.average, -0.5274E-10)
  assert approx_equal(a.standard_deviation, 0.1792E+00)
  assert a.data.origin() == (1,-4,-6)
  assert a.data.last(False) == (10,0,-3)
  assert a.data.focus() == (11,1,-2)
  assert approx_equal(a.data[:10],
    [-2.63210E-01, -4.36970E-01, -5.71930E-01, -6.09230E-01, -2.07220E-01,
     -4.15100E-01, -6.11970E-01, -7.13380E-01, -2.05500E-01, -3.60990E-01])
  assert approx_equal(a.data[40:50],
    [-4.08540E-01, -4.77320E-01, -5.16210E-01, -4.84100E-01, -2.93930E-01,
     -3.58500E-01, -4.40170E-01, -4.92110E-01, -2.19660E-01, -2.40570E-01])
  assert approx_equal(a.data[-10:],
    [-2.13550E-01, -4.87250E-01, -4.51260E-02, -2.13540E-01, -4.57070E-01,
     -6.38040E-01, -3.51570E-01, -5.98030E-01, -7.60270E-01, -7.62940E-01])
  stats = maptbx.statistics(a.data)
  assert approx_equal(stats.min(), -0.78098)
  assert approx_equal(stats.max(), 0.27233)
  assert approx_equal(stats.mean(), -0.363687)
  assert approx_equal(stats.sigma(), 0.20679)
  s = StringIO()
  a.show_summary(out=s, prefix="$")
  assert not show_diff(s.getvalue(), """\
$Title lines: 3
$   REMARKS FILENAME="cns.map"
$   REMARKS DATE:15-May-2004  02:15:56       created by user: rwgk
$   REMARKS VERSION:1.1
$Gridding:
$  n:     (24, 24, 40)
$  first: (1, -4, -6)
$  last:  (10, 0, -3)
$Total number of data points: 200
$  min:   -0.78098
$  max:   0.27233
$  mean:  -0.363687
$  sigma: 0.20679
""")
  return a

def write_xplor(map, file_name):
  iotbx.xplor.map.writer(
    file_name=file_name,
    title_lines=map.title_lines,
    unit_cell=map.unit_cell,
    gridding=map.gridding,
    data=map.data,
    average=map.average,
    standard_deviation=map.standard_deviation)

def recycle():
  for n,first,last in [[(5,3,4),(0,0,0),(3,5,6)],
                       [(4,3,5),(-1,-3,4),(6,4,5)],
                       [(3,4,5),(-2,3,0),(-2,3,0)],
                       [(3,4,5),(-2,3,0),(-2,3,3)],
                       [(3,4,5),(-2,3,0),(-2,8,0)],
                       [(3,4,5),(-2,3,0),(-2,9,0)],
                       [(3,4,5),(-2,3,0),(3,3,0)],
                       [(3,4,5),(-2,3,0),(4,3,0)]]:
    gridding = iotbx.xplor.map.gridding(
      n=n, first=first, last=last)
    flex_grid = gridding.as_flex_grid()
    data = 20000*flex.random_double(size=flex_grid.size_1d())-10000
    data.resize(flex_grid)
    stats = maptbx.statistics(data)
    iotbx.xplor.map.writer(
      file_name="tmp.map",
      title_lines=["regression test"],
      unit_cell=uctbx.unit_cell((10,20,30,80,90,100)),
      gridding=gridding,
      data=data,
      average=stats.mean(),
      standard_deviation=stats.sigma())
    read = iotbx.xplor.map.reader(file_name="tmp.map")
    assert read.title_lines == ["regression test"]
    assert read.gridding.n == n
    assert read.gridding.first == first
    assert read.gridding.last == last
    assert read.unit_cell.is_similar_to(
      uctbx.unit_cell((10,20,30,80,90,100)))
    assert eps_eq(read.average, stats.mean(), eps=1.e-4)
    assert eps_eq(read.standard_deviation, stats.sigma(), eps=1.e-4)
    assert read.data.origin() == first
    assert read.data.last(False) == last
    assert read.data.focus() == data.focus()
    assert eps_eq(read.data, data, eps=1.e-4)

def exercise_fft_map_as_xplor_map(space_group_info, n_elements=10, d_min=3):
  structure = random_structure.xray_structure(
    space_group_info,
    elements=["Si"]*n_elements,
    volume_per_atom=1000,
    min_distance=3.,
    general_positions_only=False)
  f_calc = structure.structure_factors(
    d_min=d_min, anomalous_flag=False).f_calc()
  fft_map = f_calc.fft_map()
  fft_map.as_xplor_map(
    file_name="tmp.map",
    gridding_last=[n-1 for n in fft_map.n_real()])
  read = iotbx.xplor.map.reader(file_name="tmp.map")
  assert read.title_lines == ["cctbx.miller.fft_map"]
  assert read.gridding.n == fft_map.n_real()
  assert approx_equal(flex.linear_correlation(
    read.data.as_1d(),
    fft_map.real_map_unpadded().as_1d()).coefficient(), 1)
  for first,last in [[(0,0,0),(3,5,6)],
                     [(-1,-3,4),(6,4,5)],
                     [(-2,3,0),(-2,3,0)],
                     [(-2,3,0),(-2,3,3)],
                     [(-2,3,0),(-2,8,0)],
                     [(-2,3,0),(-2,9,0)],
                     [(-2,3,0),(3,3,0)],
                     [(-2,3,0),(4,3,0)]]:
    fft_map.as_xplor_map(
      file_name="tmp.map",
      gridding_first=first,
      gridding_last=last)
    read = iotbx.xplor.map.reader(file_name="tmp.map")
    assert read.title_lines == ["cctbx.miller.fft_map"]
    assert read.gridding.n == fft_map.n_real()
    assert read.gridding.first == first
    assert read.gridding.last == last
    real_map = fft_map.real_map()
    first_p1 = [i%n for i,n in zip(first, fft_map.n_real())]
    assert eps_eq(read.data[first], real_map[first_p1], eps=1.e-4)
    last_p1 = [i%n for i,n in zip(last, fft_map.n_real())]
    assert eps_eq(read.data[last], real_map[last_p1], eps=1.e-4)
    for x in xrange(1,10):
      point = [iround(f+(l-f)*x/10.) for f,l in zip(first,last)]
      point_p1 = [i%n for i,n in zip(point, fft_map.n_real())]
      assert eps_eq(read.data[point], real_map[point_p1], eps=1.e-4)

def run():
  exercise_map_gridding()
  recycle()
  test_file_name = get_test_file_name()
  if (test_file_name is None):
    print "Skipping original CNS map test: input file not available"
  else:
    map1 = read_xplor(test_file_name)
    write_xplor(map1, "tmp.map")
    map2 = read_xplor("tmp.map")
    assert flex.max(flex.abs(map2.data-map1.data)) < 2.e-5
  exercise_fft_map_as_xplor_map(
    space_group_info=sgtbx.space_group_info("P 31"))
  print "OK"

if (__name__=="__main__"):
  run()
