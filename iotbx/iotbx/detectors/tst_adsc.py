from iotbx.detectors import adsc
import urllib

def get_test_files():
  urllib.urlretrieve('http://cci.lbl.gov/build/adsc.img','adsc.img')

def exercise_adscread():
  a = adsc.ADSCImage('adsc.img')
  a.read()
  print "Completed test read of ADSC image:"
  print "image size",a.size1,"x",a.size2,'=',a.npixels,"pixels"
  print "pixel size",a.pixel_size
  print "saturated value",a.saturation
  print "starting phi",a.osc_start

def run():
  get_test_files()
  exercise_adscread()

if __name__=="__main__":
  run()

