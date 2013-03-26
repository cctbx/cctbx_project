from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.radial_average
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import libtbx.phil
from libtbx.utils import Usage, Sorry
from cxi_user.xfel_targets import targets
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

  from iotbx.detectors.npy import NpyImage
  img = NpyImage(params.file_path)

  from spotfinder.applications.xfel import cxi_phil
  args = [params.file_path,
          "distl.detector_format_version=CXI 7.1",
          "viewer.powder_arcs.show=False",
          "viewer.powder_arcs.code=3n9c",
          ] + targets["thermolysinL785hi"]
  params_sf = cxi_phil.cxi_versioned_extract(args)
  horizons_phil = params_sf.persist.commands

  img.readHeader(horizons_phil)
  img.translate_tiles(horizons_phil)
  img.show_header()

  the_tiles = img.get_tile_manager(horizons_phil).effective_tiling_as_flex_int(
        reapply_peripheral_margin=False,encode_inactive_as_zeroes=True)
  assert len(the_tiles) == 256

  if params.beam_x is None:
    params.beam_x = img.get_beam_center_pixels_fast_slow()[0]
  if params.beam_y is None:
    params.beam_y = img.get_beam_center_pixels_fast_slow()[1]
  print "I think the beam center is (%s,%s)"%(params.beam_x, params.beam_y)

  bc = (params.beam_x,params.beam_y)

  extent = int(math.ceil(max(distance((0,0),bc),
                             distance((img.image_size_fast,0),bc),
                             distance((0,img.image_size_slow),bc),
                             distance((img.image_size_fast,img.image_size_slow),bc))))

  results = []
  for i in range(extent): results.append([])

  sys.stdout.write("Generating average...tile:")
  sys.stdout.flush()
  for tile in xrange(64):
    sys.stdout.write(" %d"%tile)
    sys.stdout.flush()

    x1,y1,x2,y2 = get_tile_coords(the_tiles,tile)
    tcx, tcy = get_tile_center(the_tiles, tile)

    for y in xrange(y1,y2):
      for x in xrange(x1,x2):
        corrected = apply_sub_pixel_metrology(tile,x,y,tcx,tcy,horizons_phil)
        if corrected is None: continue
        val = img.get_pixel_intensity((x,y))
        if val > 0:
          results[int(math.floor(distance(corrected,bc)))].append(val)

  print " Finishing..."

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

    if "%.3f"%results[i] != "nan":
      print "%.3f %.3f"%(twotheta,results[i])  #.xy format for Rex.cell.
    #print "%.3f %.3f %.3f"%(twotheta,results[i],ds[i])  # include calculated d spacings

  from pylab import scatter, show, xlim, xlabel, ylabel
  scatter(xvals,results)
  xlabel("2 theta")
  ylabel("Avg ADUs")
  show()


def get_tile_id(tiles, x, y):
    for tile in xrange(len(tiles)//4):
      x1,y1,x2,y2 = get_tile_coords(tiles, tile)
      if x <= x2 and x >= x1 and y <= y2 and y >= y1:
        return tile
    return -1

def get_tile_center(tiles, tile):
  x1,y1,x2,y2 = get_tile_coords(tiles, tile)

  cx = x1 + (distance((x2,y1),(x1,y1))/2)
  cy = y1 + (distance((x1,y2),(x1,y1))/2)
  return (cx, cy)

def distance (a,b): return math.sqrt((math.pow(b[0]-a[0],2)+math.pow(b[1]-a[1],2)))

def get_tile_coords(tiles, tile):
  """ returns x1, y1, x2, y2 """
  y1 = tiles[tile*4 + 0]
  x1 = tiles[tile*4 + 1]
  y2 = tiles[tile*4 + 2]
  x2 = tiles[tile*4 + 3]
  return (x1,y1,x2,y2)

from scitbx.matrix import col, rec, sqr, rotate_point_around_axis

def apply_sub_pixel_metrology(tile, x, y, tcx, tcy, phil):
  #return (x,y)
  if tile < 0: return None

  r = col((x, y))      # point of interest
  Ti = col((tcx,tcy))  # center of tile i

  # sub pixel translation/rotation
  ti = col((phil.integration.subpixel_joint_model.translations[tile*2],
            phil.integration.subpixel_joint_model.translations[(tile*2)+1]))
  theta = -phil.integration.subpixel_joint_model.rotations[tile] * math.pi/180

  Ri = sqr((math.cos(theta),-math.sin(theta),math.sin(theta),math.cos(theta)))

  # apply sub-pixel translation to point of interest and tile center
  rp = r + ti     # p: prime
  Tip = Ti + ti

  result = (Ri*(rp-Tip))+Tip

  return (result[0], result[1])

# Debuggind jiffy function to show a few tiles and verify metrology visually
def show_tiles(the_tiles, img, phil):
  import numpy as np
  import matplotlib.pyplot as plt

  tiles_list = (0, 1)
  #tiles_list = (2, 18)
  #tiles_list = (48, 49)
  #tiles_list = xrange(64)
  arraysx = np.array([])
  arraysy = np.array([])
  arraysz = np.array([])
  for tile in tiles_list:
    x1,y1,x2,y2 = get_tile_coords(the_tiles, tile)
    w = x2-x1
    h = y2-y1
    cx,cy = get_tile_center(the_tiles, tile)

    print "tile %d, x1 %d, y1 %d, x2 %d, y2 %d, w %d, h %d, cx %s, cy %s"%(tile, x1, y1, x2, y2, w, h, cx, cy)

    x = np.array([0]*(w*h))
    y = np.array([0]*(w*h))
    z = np.array([0]*(w*h))

    for j in xrange(h):
      for i in xrange(w):
        t_id = get_tile_id(the_tiles,x1+i,y1+j)
        if tile != t_id:
          print "bug! tile: %d, t_id %d, x %d, y %d"%(tile,t_id,x1+i,y1+j)
          return
        xt, yt = apply_sub_pixel_metrology(the_tiles,x1+i,y1+j,phil)
        x[(j*w)+i] = xt
        y[(j*w)+i] = yt
        c = img.get_pixel_intensity((j+y1,i+x1))
        z[(j*w)+i] = c

    arraysx = np.concatenate([arraysx, x])
    arraysy = np.concatenate([arraysy, y])
    arraysz = np.concatenate([arraysz, z])

  arraysz -= min(arraysz)
  arraysz /= max(arraysz)


  plt.scatter(arraysx, arraysy, c=arraysz, s=1, marker='s', cmap="gray_r", lw=0)
  plt.gca().invert_yaxis()
  #circ = plt.Circle((params.beam_x, params.beam_y), radius=200)
  #plt.gca().add_patch(circ)

  plt.show()
  return

if (__name__ == "__main__") :
  run(sys.argv[1:])
