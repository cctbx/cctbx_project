import iotbx.xplor.map
from cctbx import maptbx
from cctbx import uctbx
from cctbx.array_family import flex
from libtbx.test_utils import approx_equal
import urllib
import filecmp

def exercise_map_gridding():
  try:
    g = iotbx.xplor.map.gridding(n=(0,20,30), first=(-3,-4,-5), last=(5,4,3))
  except RuntimeError, e:
    assert str(e) == "Illegal xplor map gridding for dimension X: " \
                   + "gridding=0, first=-3, last=5"
  g = iotbx.xplor.map.gridding(n=(10,20,30), first=(-3,-4,-5), last=(5,4,3))
  fg = g.as_flex_grid()
  assert fg.origin() == g.first
  assert fg.last(0) == g.last
  assert fg.last(1) == (6,5,4)
  assert fg.all() == (9,9,9)
  assert not fg.is_padded()

def get_test_files():
  urllib.urlretrieve(
    'http://cci.lbl.gov/cctbx_downloads/regression/iotbx/NSFN_C2221.xplor',
    'NSFN_C2221.xplor')

def read_xplor(file_name):
  a = iotbx.xplor.map.reader(file_name=file_name)
  assert a.title_lines == [' REMARKS FILENAME=""',
                           ' REMARKS Phenix Xarray to CNS map format']
  assert a.gridding.n == (24,120,54)
  assert a.gridding.first == (0,0,0)
  assert a.gridding.last == (24,120,54)
  assert a.unit_cell.is_similar_to(
    uctbx.unit_cell((32.042, 175.362, 79.663, 90, 90, 90)))
  assert approx_equal(a.average, 0)
  assert approx_equal(a.standard_deviation, 1)
  d = a.data
  assert approx_equal(d[0:5], (-0.284546,-0.60108,-1.11654,-1.17415,-0.827963))
  assert a.data.origin() == (0,0,0)
  assert a.data.last(0) == (24,120,54)
  assert a.data.last() == (25,121,55)
  assert a.data.focus() == (25,121,55)
  assert approx_equal(flex.min(a.data), -1.86868)
  assert approx_equal(flex.max(a.data), 4.18869)
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
    data = flex.random_double(size=flex_grid.size_1d())
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
    assert approx_equal(read.average, stats.mean(), eps=2.e-4)
    assert approx_equal(read.standard_deviation, stats.sigma(), eps=2.e-4)
    assert read.data.origin() == first
    assert read.data.last(0) == last
    assert read.data.focus() == data.focus()
    assert flex.max(flex.abs(read.data-data)) < 2.e-5

def run():
  exercise_map_gridding()
  recycle()
  get_test_files()
  map1 = read_xplor("NSFN_C2221.xplor")
  write_xplor(map1, "comparison.xplor")
  assert filecmp.cmp("NSFN_C2221.xplor", "comparison.xplor")
  map2 = read_xplor("comparison.xplor")
  assert flex.max(flex.abs(map2.data-map1.data)) < 2.e-5
  print "OK"

if (__name__=="__main__"):
  run()
