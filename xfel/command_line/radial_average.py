from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.radial_average
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import libtbx.phil
from libtbx.utils import Usage, Sorry
import sys
import os
import math
import numpy as np

master_phil = libtbx.phil.parse("""
  file_path = None
    .type = str
  beam_x = None
    .type = float
  beam_y = None
    .type = float
""")

def run (args) :
  user_phil = []
  # TODO: replace this stuff with iotbx.phil.process_command_line_with_files
  # as soon as I can safely modify it
  for arg in args :
    if (not "=" in arg) :
      try :
        user_phil.append(libtbx.phil.parse("""file_path=%s""" % arg))
      except ValueError, e :
        raise Sorry("Unrecognized argument '%s'" % arg)
    else :
      try :
        user_phil.append(libtbx.phil.parse(arg))
      except RuntimeError, e :
        raise Sorry("Unrecognized argument '%s' (error: %s)" % (arg, str(e)))
  params = master_phil.fetch(sources=user_phil).extract()
  if params.file_path is None or not os.path.isfile(params.file_path) :
    master_phil.show()
    raise Usage("file_path must be defined (either file_path=XXX, or the path alone).")

  from rstbx.command_line.slip_viewer import ImageFactory
  try:
    img = ImageFactory(params.file_path)
  except Exception,e:
    raise Sorry("Couldn't load the image:" + e)



  if params.beam_x is None:
    params.beam_x = img.get_beam_center_pixels_fast_slow()[0]
  if params.beam_y is None:
    params.beam_y = img.get_beam_center_pixels_fast_slow()[1]

  print "I think the beam center is (%s,%s)"%(params.beam_x, params.beam_y)
  def distance (a,b): return math.sqrt((math.pow(b[0]-a[0],2)+math.pow(b[1]-a[1],2)))

  bc = (params.beam_x,params.beam_y)
  extent = int(math.ceil(max(distance((0,0),bc),
                             distance((img.image_size_fast,0),bc),
                             distance((0,img.image_size_slow),bc),
                             distance((img.image_size_fast,img.image_size_slow),bc))))

  results = []
  for i in range(extent): results.append([])

  sys.stdout.write("Generating average...")
  sys.stdout.flush()
  window = int(img.image_size_slow/10)
  for y in range(img.image_size_slow):
    for x in range(img.image_size_fast):
      val = img.get_pixel_intensity((x,y))
      if val > 0:
        results[int(math.floor(distance((x,y),bc)))].append(val)

    if(y%window==0):
      sys.stdout.write("%d%%..."%(100*y/(window*10)))
      sys.stdout.flush()
      if y != 0: pass

  print "Finishing..."

  xvals = np.ndarray((extent,),float)
  #ds = np.ndarray((extent,),float)
  for i in range(extent):
    results[i] = np.mean(results[i])
    d_in_mm = i * img.pixel_resolution
    twotheta = math.atan(d_in_mm/img.distance)*180/math.pi
    xvals[i] = twotheta

    #if twotheta==0:
      #ds[i] = 1000
    #else:
      #ds[i] = img.wavelength/(math.sin(math.pi*twotheta/180))

    print "%.3f %.3f"%(twotheta,results[i])  #.xy format for Rex.cell.  TODO remove nan values
    #print "%.3f %.3f %.3f"%(twotheta,results[i],ds[i])  # include calculated d spacings

  from pylab import scatter, show, xlim, xlabel, ylabel
  scatter(xvals,results)
  xlabel("2 theta")
  ylabel("Avg ADUs")
  show()

if (__name__ == "__main__") :
  run(sys.argv[1:])
