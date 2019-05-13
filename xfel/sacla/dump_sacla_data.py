from __future__ import division
from __future__ import print_function
import sys, os, dxtbx
from xfel.cxi.cspad_ana import cspad_tbx
from libtbx import easy_pickle
from libtbx import easy_mp
from xfel.command_line.cxi_image2pickle import crop_image_pickle

# Jiffy script to dump SACLA data processed by Cheetah into image pickles.  Usage:
# libtbx.python dump_sacla_data.py <path to h5 file> <destination directory for pickles>.
# Uses 4 processors, hardcoded at the end of the file

data_path = sys.argv[1]
dest_dir = sys.argv[2]

data = dxtbx.load(data_path)
detector = data.get_detector()
distance = detector[0].get_directed_distance()
beam = data.get_beam()
wavelength = beam.get_wavelength()
pixel_size = detector[0].get_pixel_size()[0]
beam_x, beam_y = detector[0].get_beam_centre_px(beam.get_s0())
beam_x *= pixel_size
beam_y += 0 #earlier work required a 3 pixel shift based on powder pattern/fit to unit cell
beam_y *= pixel_size
overload = detector[0].get_trusted_range()[1]
from xfel.cxi.cspad_ana.cspad_tbx import xpp_active_areas
active_areas = xpp_active_areas["Sacla.MPCCD.8tile"]["active_areas"]
# the active areas are already determined for the cropped size
# (ran once without active areas, then measured cropped active areas on image viewer)

dest_base = os.path.basename(os.path.splitext(data_path)[0])

def do_work(img_no):
  n_fails = 0
  while True:
    try:
      raw_data = data.get_raw_data(img_no)
      break
    except (KeyError, ValueError):
      n_fails +=1
      print("Fail to read, attempt number", n_fails)
      if n_fails > 100:
        raise Exception("Couldn't read the data")
    import time; time.sleep(n_fails * 0.1)

  imgdict = cspad_tbx.dpack(data=raw_data,
     distance=distance,
     pixel_size=pixel_size,
     wavelength=wavelength,
     beam_center_x=beam_x,
     beam_center_y=beam_y,
     ccd_image_saturation=overload,
     saturated_value=overload,
     address="Sacla.MPCCD.8tile",
     active_areas=active_areas
     )
  imgdict = crop_image_pickle(imgdict,
    preserve_active_areas_even_though_cropping_would_invalidate_them=True)

  dest_path = os.path.join(dest_dir, dest_base + "_%06d.pickle"%img_no)
  print("Saving image", img_no, "to", dest_path)
  easy_pickle.dump(dest_path, imgdict)

easy_mp.pool_map(
            args=list(range(data.get_num_images())),
            func=do_work,
            processes=4)
