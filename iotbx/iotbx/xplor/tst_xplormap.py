from iotbx.xplor import map_gridding, XplorMap
from cctbx import uctbx
from libtbx.test_utils import approx_equal
import urllib
import filecmp

def exercise_map_gridding():
  try:
    g = map_gridding(n=(0,20,30), first=(-3,-4,-5), last=(5,4,3))
  except RuntimeError, e:
    assert str(e) == "Illegal xplor map gridding for dimension X: " \
                   + "gridding=0, first=-3, last=5"
  g = map_gridding(n=(10,20,30), first=(-3,-4,-5), last=(5,4,3))
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
  a = XplorMap().read(file_name=file_name)
  assert a.title == [' REMARKS FILENAME=""',
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
  return a

def write_xplor(map, file_name):
  a = XplorMap()
  a.title = map.title
  a.gridding = map.gridding
  a.unit_cell = map.unit_cell
  a.data = map.data
  a.avarage = map.average
  a.standard_deviation = map.standard_deviation
  a.write(file_name=file_name)

def run():
  exercise_map_gridding()
  get_test_files()
  map = read_xplor("NSFN_C2221.xplor")
  write_xplor(map, "comparison.xplor")
  assert filecmp.cmp("NSFN_C2221.xplor", "comparison.xplor")
  print "OK"

if (__name__=="__main__"):
  run()
