from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 05/19/2015
Description : Converts raw image to pickle files; crops or pads pickle to place
              beam center into center of image; masks out beam stop. Adapted in
              part from cxi_image2pickle.py by Aaron Brewster.
'''

import os
import sys
from cStringIO import StringIO
import math

from scitbx.array_family import flex

import dxtbx
from cPickle import load
from libtbx import easy_pickle as ep
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp
from libtbx.easy_mp import parallel_map

import iota_cmd as cmd

class Capturing(list):
  """ Class used to capture stdout from cctbx.xfel objects. Saves output in
  appendable list for potential logging.
  """
  def __enter__(self):
    self._stdout = sys.stdout
    self._stderr = sys.stderr
    sys.stdout = self._stringio_stdout = StringIO()
    sys.stderr = self._stringio_stderr = StringIO()
    return self
  def __exit__(self, *args):
    self.extend(self._stringio_stdout.getvalue().splitlines())
    sys.stdout = self._stdout
    self.extend(self._stringio_stderr.getvalue().splitlines())
    sys.stderr = self._stderr

def load_image(img):
  """ Reads raw image file (tested for MCCD so far, only) and extracts data for
      conversion into pickle format.

      input: img - raw image file (abs. path)

      output: data - image file info (w/ raw data) for pickling
  """

  try:
    with Capturing() as junk_output:
      loaded_img = dxtbx.load(img)
  except IOError:
    loaded_img = None
    pass

  if loaded_img is not None:
    raw_data   = loaded_img.get_raw_data()
    detector   = loaded_img.get_detector()[0]
    beam       = loaded_img.get_beam()
    scan       = loaded_img.get_scan()

    distance   = detector.get_distance()
    pixel_size = detector.get_pixel_size()[0]
    overload   = detector.get_trusted_range()[1]
    wavelength = beam.get_wavelength()
    beam_x     = detector.get_beam_centre(beam.get_s0())[0]
    beam_y     = detector.get_beam_centre(beam.get_s0())[1]

    if scan is None:
      timestamp = None
    else:
      msec, sec = math.modf(scan.get_epochs()[0])
      timestamp = evt_timestamp((sec,msec))

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
  else:
    data = None

  return data


def check_image(img):

  try:
    with Capturing() as suppressed_junk:
      loaded_img = dxtbx.load(img)
  except IOError:
    loaded_img = None
    pass

  if loaded_img is not None:
    scan = loaded_img.get_scan()
    if scan is None:
      detector   = loaded_img.get_detector()[0]
      beam       = loaded_img.get_beam()
      beam_x     = detector.get_beam_centre(beam.get_s0())[0]
      beam_y     = detector.get_beam_centre(beam.get_s0())[1]
      if int(beam_x) == int(beam_y):
        verdict = 'converted pickle'
      else:
        verdict = 'raw pickle'
    else:
      verdict = 'image'
  else:
    print img
    verdict = None

  return verdict


def square_pickle(data, square='crop'):
  """ A function to crop the image pickle to a square such that the beam center
      is in the center of image (mostly copied from cxi_image2pickle.py)

      input: data - image data

      output: data - amended image data
  """

  # only one active area is allowed, and it should be the size of the image.

  test = flex.int([0,0,data['SIZE1'],data['SIZE2']]) == data['ACTIVE_AREAS']
  assert test.all_eq(True)

  # retrieve parameters from the dictionary
  pixel_size = data['PIXEL_SIZE']
  beam_x = int(round(data['BEAM_CENTER_X']/pixel_size))
  beam_y = int(round(data['BEAM_CENTER_Y']/pixel_size))
  width = data['SIZE1']
  height = data['SIZE2']
  pixels = data['DATA']
  right = beam_x
  left = width - beam_x
  top = beam_y
  bottom = height - beam_y
  new_pixels = pixels.deep_copy()

  # the new image will be twice the size as the smallest distance between the
  # beam center and one of the edges of the images
  if square == 'crop':
    new_half_size = min([right, left, top, bottom])
    new_size = new_half_size * 2
    new_size_x = new_size_y = new_size
    min_x = beam_x - new_half_size
    min_y = beam_y - new_half_size
    new_beam_x = data['BEAM_CENTER_X'] - min_x * pixel_size
    new_beam_y = data['BEAM_CENTER_Y'] - min_y * pixel_size
    new_pixels = pixels[min_y:min_y+new_size,min_x:min_x+new_size]
    assert new_pixels.focus()[0] == new_pixels.focus()[1]

  # the new image will be twice the size as the LARGEST distance between the
  # beam center and one of the edges of the images
  elif square == 'pad':
    new_half_size = max([right, left, top, bottom])
    new_size = new_half_size * 2
    new_size_x = new_size_y = new_size
    min_x = 0
    min_y = 0
    delta_x = new_half_size - beam_x
    delta_y = new_half_size - beam_y
    new_beam_x = data['BEAM_CENTER_X'] + delta_x * pixel_size
    new_beam_y = data['BEAM_CENTER_Y'] + delta_y * pixel_size

    new_pixels.resize(flex.grid(new_size, new_size))
    new_pixels.fill(-2)

    for y in xrange(pixels.focus()[0]):
      for x in xrange(pixels.focus()[1]):
        new_pixels[y + delta_y, x + delta_x] = pixels[y, x]
    assert new_pixels.focus()[0] == new_pixels.focus()[1]
#   else:
#     new_size_x = width
#     new_size_y = height

  # save the results
  data['DATA'] = new_pixels
  data['SIZE1'] = new_size_x
  data['SIZE2'] = new_size_y
  data['BEAM_CENTER_X'] = new_beam_x
  data['BEAM_CENTER_Y'] = new_beam_y
  data['ACTIVE_AREAS'] = flex.int([0 , 0, new_size, new_size])

  return data


def create_filenames(img_file, dest_dir):
  """ Generates filenames for converted files

      input: img_file - original image file (absolute path)
             dest_dir - destination folder (absolute path)

      output: masked_file - filename (w/ abs. path) for converted image file
  """


  filename = os.path.basename(img_file)
  filename_no_ext = filename.split('.')[0]

  masked_filename = filename_no_ext + "_prep.pickle"
  masked_file = os.path.join(dest_dir, masked_filename)

  return masked_file

def mask_image(data):
  """ Identifies beamstop shadow and sets pixels inside to -2 (to be ignored by
      processing software). Clunky now (merely finds all pixels with intensity
      less than 0.4 * image average intensity), to be refined later.

      input: data - image data

      output: data - modified image data
  """

  img_raw_bytes = data['DATA']

  img_thresh = int((0.45*flex.mean(img_raw_bytes.as_double())))
  beam_stop_sel = img_raw_bytes <= img_thresh
  img_masked = img_raw_bytes.set_selected(beam_stop_sel, -2)

  data['DATA'] = img_masked
  return data


def convert_image(img_in, img_out, square, beamstop=True):
  """ Converts images into pickle format; crops and masks out beamstop if
      selected

      input: img - image filename (absolute path)
             dest_dir - destination folder for converted image
             square - toggle making image a square w/ beam center in center of
                      image. Only one option now (to crop); padding may or may
                      not be added later.
             beamstop - toggle masking out beamstop shadow

      output: writes down a file in pickle format, preserving filename
  """

  img_data = load_image(img_in)

  if square != "None":
    img_data = square_pickle(img_data, square)

  if beamstop:
    img_data = mask_image(img_data)

  if not os.path.isdir(os.path.dirname(img_out)):
    os.makedirs(os.path.dirname(img_out))
  ep.dump(img_out, img_data)


def run_wrapper(img_in):
  """ Multiprocessor wrapper, for testing module
  """

  img_out = create_filenames(img_in, dest_folder)
  prog_count = img_files.index(img_in)
  n_img = len(img_files)

  gs_prog = cmd.ProgressBar(title='CONVERTING IMAGES')
  if prog_count < n_img:
    prog_step = 100 / n_img
    gs_prog.update(prog_count * prog_step, prog_count)
  else:
    gs_prog.finished()

  return convert_image(img_in, img_out, square, beamstop)


# ============================================================================ #

if __name__ == "__main__":

  data_folder = os.path.abspath(sys.argv[1])
  dest_folder = os.path.abspath(sys.argv[2])
  img_files = [os.path.join(data_folder, f) for f in os.listdir(data_folder) if os.path.isfile(os.path.join(data_folder, f))]

  beamstop = True

  if sys.argv[3] == '-c':
    square = 'crop'
  elif sys.argv[3] == '-p':
    square = 'pad'
  elif sys.argv[3] == '-x':
    square = "None"
    beamstop = False
  else:
    square = "None"
    if sys.argv[3] == '-t':
      print check_image(img_files[0])
      sys.exit()



  cmd.Command.start("Converting {} images".format(len(img_files)))
  parallel_map(iterable  = img_files,
               func      = run_wrapper,
               processes = 12)
  cmd.Command.end("Converting {} images -- DONE ".format(len(img_files)))
