from __future__ import division
# LIBTBX_SET_DISPATCHER_NAME cxi.image2pickle

# Convert images of any extant format to pickle files suitable for processing with
# cxi.index.  Note, oscillation values are not preserved.

import dxtbx, sys, os, math
import libtbx.option_parser
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp
from libtbx import easy_pickle
from libtbx.utils import Usage

def run(argv=None):
  if argv is None:
    argv = sys.argv[1:]

  command_line = (libtbx.option_parser.option_parser(
    usage="%s [-v] " % libtbx.env.dispatcher_name)
                  .option(None, "--verbose", "-v",
                          action="store_true",
                          default=False,
                          dest="verbose",
                          help="Print more information about progress")
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
    destpath = os.path.join(os.path.dirname(imgpath), os.path.splitext(os.path.basename(imgpath))[0] + ".pickle")

    img = dxtbx.load(imgpath)

    detector = img.get_detector()
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


    beam = img.get_beam()
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
        raise Usage("Can't get beam x position from image. Override with -x")
      if command_line.options.beam_center_y is None:
        raise Usage("Can't get beam y position from image. Override with -y")
      beam_x = command_line.options.beam_center_x * pixel_size
      beam_y = command_line.options.beam_center_y * pixel_size
    else:
      if command_line.options.beam_center_x is None:
        beam_x = detector.get_beam_centre(beam.get_s0())[0]
      else:
        beam_x = command_line.options.beam_center_x

      if command_line.options.beam_center_y is None:
        beam_y = detector.get_beam_centre(beam.get_s0())[1]
      else:
        beam_x = command_line.options.beam_center_x

    scan = img.get_scan()
    if scan is None:
      timestamp = None
    else:
      msec, sec = math.modf(scan.get_epochs()[0])
      timestamp = evt_timestamp((sec,msec))

    data = dpack(data=img.get_raw_data(),
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

    easy_pickle.dump(destpath, data)

if (__name__ == "__main__") :
  run(sys.argv[1:])
