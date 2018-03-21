from __future__ import absolute_import, division
from __future__ import print_function
# LIBTBX_SET_DISPATCHER_NAME dxtbx.image2pickle
# LIBTBX_SET_DISPATCHER_NAME cxi.image2pickle

# Convert images of any extant format to pickle files suitable for processing with
# cxi.index.  Note, oscillation values are not preserved.

import dxtbx, sys, os, math
import libtbx.option_parser
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp
from libtbx import easy_pickle
from libtbx.utils import Usage

def crop_image_pickle(data,preserve_active_areas_even_though_cropping_would_invalidate_them=False):
  """
  Given an image pickle dictionary, crop the pixels such that the beam center is as close
  as possile to the image center.  Then adjust SIZE1/SIZE2, ACTIVE_AREAS and the beam
  center accordingly.
  @param data The image dictionary of interest
  """
  # only one active area is allowed, and it should be the size of the image.
  from scitbx.array_family import flex
  if preserve_active_areas_even_though_cropping_would_invalidate_them is False:
    test = flex.int([0,0,data['SIZE1'],data['SIZE2']]) == data['ACTIVE_AREAS']
    assert test.all_eq(True)

  # retrieve parameters from the dictionary
  pixel_size = data['PIXEL_SIZE']
  beam_x = int(round(data['BEAM_CENTER_X']/pixel_size))
  beam_y = int(round(data['BEAM_CENTER_Y']/pixel_size))
  width = data['SIZE1']
  height = data['SIZE2']
  pixels = data['DATA']

  # the new image will be twice the size as the smallest distance between the
  # beam center and one of the edges of the images
  new_half_size = min([beam_x, width-beam_x, beam_y, height-beam_y])
  new_size = new_half_size * 2
  min_x = beam_x-new_half_size
  min_y = beam_y-new_half_size
  pixels = pixels[min_y:min_y+new_size,min_x:min_x+new_size]
  assert pixels.focus()[0] == pixels.focus()[1]

  # save the results
  data['DATA'] = pixels
  data['SIZE1'] = new_size
  data['SIZE2'] = new_size
  data['BEAM_CENTER_X'] -= min_x * pixel_size
  data['BEAM_CENTER_Y'] -= min_y * pixel_size
  if preserve_active_areas_even_though_cropping_would_invalidate_them is False:
    data['ACTIVE_AREAS'] = flex.int([0,0,new_size,new_size])

  return data

def run(argv=None):
  if argv is None:
    argv = sys.argv[1:]

  command_line = (libtbx.option_parser.option_parser(
    usage="%s [-v] [-c] [-s] [-w wavelength] [-d distance] [-p pixel_size] [-x beam_x] [-y beam_y] [-o overload] files" % libtbx.env.dispatcher_name)
                  .option(None, "--verbose", "-v",
                          action="store_true",
                          default=False,
                          dest="verbose",
                          help="Print more information about progress")
                  .option(None, "--crop", "-c",
                          action="store_true",
                          default=False,
                          dest="crop",
                          help="Crop the image such that the beam center is in the middle")
                  .option(None, "--skip_converted", "-s",
                          action="store_true",
                          default=False,
                          dest="skip_converted",
                          help="Skip converting if an image already exist that matches the destination file name")
                  .option(None, "--wavelength", "-w",
                          type="float",
                          default=None,
                          dest="wavelength",
                          help="Override the image's wavelength (angstroms)")
                  .option(None, "--distance", "-d",
                          type="float",
                          default=None,
                          dest="distance",
                          help="Override the detector distance (mm)")
                  .option(None, "--pixel_size", "-p",
                          type="float",
                          default=None,
                          dest="pixel_size",
                          help="Override the detector pixel size (mm)")
                  .option(None, "--beam_x", "-x",
                          type="float",
                          default=None,
                          dest="beam_center_x",
                          help="Override the beam x position (pixels)")
                  .option(None, "--beam_y", "-y",
                          type="float",
                          default=None,
                          dest="beam_center_y",
                          help="Override the beam y position (pixels)")
                  .option(None, "--overload", "-o",
                          type="float",
                          default=None,
                          dest="overload",
                          help="Override the detector overload value (ADU)")
                  ).process(args=argv)

  paths = command_line.args
  if len(paths) <= 0:
    raise Usage("No files specified")

  for imgpath in paths:
    if command_line.options.verbose:
      print("Reading %s"%(imgpath))

    try:
      img = dxtbx.load(imgpath)
    except IOError:
      img = None
      pass

    if img is None:
      import numpy as np
      try:
        raw_data = np.loadtxt(imgpath)

        from scitbx.array_family import flex
        raw_data = flex.double(raw_data.astype(np.double))
      except ValueError:
        raise Usage("Couldn't load %s, no supported readers"%imgpath)

      detector = None
      beam = None
      scan = None
      is_multi_image = False
    else:
      try:
        raw_data = img.get_raw_data()
        is_multi_image = False
      except TypeError:
        raw_data = img.get_raw_data(0)
        is_multi_image = True
      detector = img.get_detector()
      beam = img.get_beam()
      scan = img.get_scan()


    if detector is None:
      if command_line.options.distance is None:
        raise Usage("Can't get distance from image. Override with -d")
      if command_line.options.pixel_size is None:
        raise Usage("Can't get pixel size from image. Override with -p")
      if command_line.options.overload is None:
        raise Usage("Can't get overload value from image. Override with -o")
      distance = command_line.options.distance
      pixel_size = command_line.options.pixel_size
      overload = command_line.options.overload
    else:
      detector = detector[0]
      if command_line.options.distance is None:
        distance   = detector.get_distance()
      else:
        distance = command_line.options.distance

      if command_line.options.pixel_size is None:
        pixel_size = detector.get_pixel_size()[0]
      else:
        pixel_size = command_line.options.pixel_size

      if command_line.options.overload is None:
        overload   = detector.get_trusted_range()[1]
      else:
        overload = command_line.options.overload

    if beam is None:
      if command_line.options.wavelength is None:
        raise Usage("Can't get wavelength from image. Override with -w")
      wavelength = command_line.options.wavelength
    else:
      if command_line.options.wavelength is None:
        wavelength = beam.get_wavelength()
      else:
        wavelength = command_line.options.wavelength

    if beam is None and detector is None:
      if command_line.options.beam_center_x is None:
        print("Can't get beam x position from image. Using image center. Override with -x")
        beam_x = raw_data.focus()[0] * pixel_size
      else:
        beam_x = command_line.options.beam_center_x * pixel_size

      if command_line.options.beam_center_y is None:
        print("Can't get beam y position from image. Using image center. Override with -y")
        beam_y = raw_data.focus()[1] * pixel_size
      else:
        beam_y = command_line.options.beam_center_y * pixel_size
    else:
      if command_line.options.beam_center_x is None:
        beam_x = detector.get_beam_centre(beam.get_s0())[0]
      else:
        beam_x = command_line.options.beam_center_x * pixel_size

      if command_line.options.beam_center_y is None:
        beam_y = detector.get_beam_centre(beam.get_s0())[1]
      else:
        beam_y = command_line.options.beam_center_y * pixel_size

    if scan is None:
      timestamp = None
    else:
      msec, sec = math.modf(scan.get_epochs()[0])
      timestamp = evt_timestamp((sec,msec))

    if is_multi_image:
      for i in xrange(img.get_num_images()):
        save_image(command_line, imgpath, scan, img.get_raw_data(i), distance, pixel_size, wavelength, beam_x, beam_y, overload, timestamp, image_number = i)
    else:
      save_image(command_line, imgpath, scan, raw_data, distance, pixel_size, wavelength, beam_x, beam_y, overload, timestamp)

def save_image(command_line, imgpath, scan, raw_data, distance, pixel_size, wavelength, beam_x, beam_y, overload, timestamp, image_number = None):
  if image_number is None:
    destpath = os.path.join(os.path.dirname(imgpath), os.path.splitext(os.path.basename(imgpath))[0] + ".pickle")
  else:
    destpath = os.path.join(os.path.dirname(imgpath), os.path.splitext(os.path.basename(imgpath))[0] + "%05d.pickle"%image_number)
  if command_line.options.skip_converted and os.path.isfile(destpath):
    if command_line.options.verbose:
      print("Skipping %s, file exists"%imgpath)
      return

  data = dpack(data=raw_data,
               distance=distance,
               pixel_size=pixel_size,
               wavelength=wavelength,
               beam_center_x=beam_x,
               beam_center_y=beam_y,
               ccd_image_saturation=overload,
               saturated_value=overload,
               timestamp=timestamp
               )

  if scan is not None:
    osc_start, osc_range = scan.get_oscillation()
    if osc_start != osc_range:
      data['OSC_START'] = osc_start
      data['OSC_RANGE'] = osc_range

      data['TIME'] = scan.get_exposure_times()[0]

  if command_line.options.crop:
    data = crop_image_pickle(data)

  if command_line.options.verbose:
    print("Writing", destpath)

  easy_pickle.dump(destpath, data)

if (__name__ == "__main__") :
  run(sys.argv[1:])
