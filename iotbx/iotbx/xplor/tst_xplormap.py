from iotbx.xplor import XplorMap
import urllib

def get_test_files():
  urllib.urlretrieve('http://cci.lbl.gov/build/NSFN_C2221.xplor','NSFN_C2221.xplor')

def read_xplor(f):
  print "Test of Xplor map read:"
  a = XplorMap().read(f)
  print a.title
  print a.sections
  print ["%.3f"%v for v in a.unitcell.parameters()]
  print a.order
  print a.average
  print a.stddev
  d = a.data
  print tuple(d[0:5])
  print a.data.focus()
  return a
  
def write_xplor(map,f):
  print "Test of Xplor map write:"
  a = XplorMap()
  a.title=map.title
  a.sections=map.sections
  a.unitcell=map.unitcell
  a.order=map.order
  a.data=map.data
  a.write(f)
  
  
def run():
  #get_test_files()
  map = read_xplor('NSFN_C2221.xplor')
  write_xplor(map,'comparison.xplor')

if __name__=="__main__":
  run()

