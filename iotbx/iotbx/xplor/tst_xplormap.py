from iotbx.xplor import XplorMap
from libtbx.test_utils import approx_equal
import urllib
import filecmp

def get_test_files():
  urllib.urlretrieve(
    'http://cci.lbl.gov/cctbx_downloads/regression/iotbx/NSFN_C2221.xplor',
    'NSFN_C2221.xplor')

def read_xplor(f):
  a = XplorMap().read(f)
  assert a.title == [' REMARKS FILENAME=""',
                     ' REMARKS Phenix Xarray to CNS map format']
  assert a.sections == [24, 0, 24, 120, 0, 120, 54, 0, 54]
  assert ["%.3f"%v for v in a.unitcell.parameters()] \
      == ['32.042', '175.362', '79.663', '90.000', '90.000', '90.000']
  assert a.order == "ZYX"
  assert approx_equal(a.average, 0)
  assert approx_equal(a.stddev, 1)
  d = a.data
  assert approx_equal(d[0:5], (-0.284546,-0.60108,-1.11654,-1.17415,-0.827963))
  assert a.data.focus() == (25, 121, 55)
  return a

def write_xplor(map,f):
  a = XplorMap()
  a.title=map.title
  a.sections=map.sections
  a.unitcell=map.unitcell
  a.order=map.order
  a.data=map.data
  a.write(f)

def run():
  get_test_files()
  map = read_xplor('NSFN_C2221.xplor')
  write_xplor(map,'comparison.xplor')
  assert filecmp.cmp('NSFN_C2221.xplor','comparison.xplor')
  print "OK"

if __name__=="__main__":
  run()
