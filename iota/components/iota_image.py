from __future__ import division, print_function, absolute_import
from past.builtins import range

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 11/05/2018
Description : Creates image object. If necessary, converts raw image to pickle
              files; crops or pads pickle to place beam center into center of
              image; masks out beam stop. (Adapted in part from
              cxi_image2pickle.py by Aaron Brewster.) If selected, estimates gain
              (from dials.estimate_gain by Richard Gildea). If selected, checks image
              for diffraction. Creates and runs instances of integrator (and selector,
              in case of cctbx.xfel) objects.
'''

import os
import math
import numpy as np

try:  # for Py3 compatibility
    import itertools.ifilter as filter
except ImportError:
    pass


from scitbx.array_family import flex

import dxtbx
from libtbx import easy_pickle as ep
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp

import iota.components.iota_utils as util
from iota.components.iota_base import SingleImageBase, ImageImporterBase


class OldSingleImage(SingleImageBase):
  ''' SingleImage object for use with the old HA14 processing (with
  conversion and all those goodies) '''
  def __init__(self, imgpath, idx=None):
    SingleImageBase.__init__(self, imgpath=imgpath, idx=idx)
    self.raw_img = imgpath
    self.conv_img = None
    self.gain = 1.0

    # Grid search parameters
    self.grid = []
    self.gs_results = []

    # Append DISTL spot-finding parameters to final info object
    self.final['sih'] = 0       # Signal height
    self.final['sph'] = 0       # Spot height
    self.final['spa'] = 0       # Spot area

  def generate_grid(self, gs_type='brute_force', sig_search=False,
                    hr=0, ar=0, hmed=0, amed=0):
    ''' HA14-specific function to run image triage and generate grid search '''

    if gs_type == 'smart':
      hrange = 1
      arange = 1
    elif gs_type in ('none', 'no_grid_search'):
      hrange = 0
      arange = 0
    else:
      hrange = hr
      arange = ar
    grid_points = []

    h_min = hmed - hrange
    h_max = hmed + hrange
    a_min = amed - arange
    a_max = amed + arange
    h_std = hr
    a_std = ar

    for spot_area in range(a_min, a_max + 1):
      for spot_height in range(h_min, h_max + 1):
        if sig_search:
          if spot_height >= 1 + h_std:
            sigs = range(spot_height - h_std, spot_height + 1)
          elif spot_height < 1 + h_std:
            sigs = range(1, spot_height + 1)
          elif spot_height == 1:
            sigs = [1]
        else:
          sigs = [spot_height]

        for sig_height in sigs:
          if (spot_area, spot_height, sig_height) not in grid_points:
            grid_points.append((spot_area, spot_height, sig_height))
            self.grid.append(
              {'sih': sig_height,
               'sph': spot_height,
               'spa': spot_area,
               'a': 0, 'b': 0, 'c': 0,
               'alpha': 0, 'beta': 0, 'gamma': 0,
               'sg': '', 'strong': 0, 'res': 0, 'lres': 0,
               'mos': 0, 'epv': 0, 'info': '', 'ok': True}
            )

class OldImageImporter(ImageImporterBase):
  ''' Image importer for old-skool HA14 processing (outputs an image pickle) '''

  def __init__(self, init):
    ImageImporterBase.__init__(self, init=init)
    self.params = init.params
    self.modify = True
    self.conv_base = util.set_base_dir('converted_pickles',
                                       out_dir=self.params.output)

  def instantiate_image_object(self, filepath, idx=None):
    ''' Instantiate a SingleImage object for current backend
    :param filepath: path to image file
    '''
    self.img_object = OldSingleImage(imgpath=filepath, idx=idx)

  def image_triage(self):
    # Triage image (i.e. check for usable diffraction, using selected method)
    from iota.components.iota_cctbx_ha14 import Triage
    triage = Triage(self.img_object.conv_img, self.params)
    fail, log_entry, hmed, amed = triage.triage_image()

    if not fail:
      self.img_object.status = 'triaged'
    else:
      self.img_object.status = fail

    self.img_object.fail = fail
    self.img_object.log_info.append(log_entry)

    return hmed, amed

  def load_image_file(self, filepath):
    """ Reads a single image file (e.g. a CBF or MCCD) and returns raw data,
    experiment information, and an image type (e.g. pickle, raw image, or none)
    :param img: path to image file
    :return: data: raw image data and experiment info
             img_type: which type of image this was, or None if not loaded
    """

    error = None
    try:
      with util.Capturing() as junk_output:
        loaded_img = dxtbx.load(filepath)
    except Exception as e:
      error = 'IOTA IMPORTER ERROR: DXTBX failed! {}'.format(e)
      print (e)
      loaded_img = None

    # Extract image information
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

      # Assemble datapack
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
          data['OSC_START'] = 0 #osc_start
          data['OSC_RANGE'] = 0 #osc_start
          data['TIME'] = scan.get_exposure_times()[0]

      # Populate image object information
      self.img_object.final['pixel_size'] = pixel_size
      self.img_object.final['img_size'] = (data['SIZE1'], data['SIZE2'])
      self.img_object.final['beamX'] = beam_x
      self.img_object.final['beamY'] = beam_y
      self.img_object.final['gain'] = detector.get_gain()
      self.img_object.final['distance'] = distance
      self.img_object.final['wavelength'] = wavelength

    else:
      data = None

    return data, error

  def apply_mask_from_file(self, data, mask_file, invert=False):
    ''' Read a mask from file (should be a bool array) and set every pixel
    under mask to -2 (i.e. will be ignored by processing)

    :param data: raw image data
    :param mask_file: filepath to the mask
    :param invert: invert boolean array from True to False (need for HA14)
    :return: masked data '''

    img_raw_bytes = data['DATA']
    error = None

    try:
      full_mask = ep.load(mask_file)
      if type(full_mask) == tuple:
        full_mask = full_mask[0]

      if invert:
        mask = full_mask == False
      else:
        mask = full_mask
      if mask is not None:
        img_masked = img_raw_bytes.set_selected(mask, -2)
        data['DATA'] = img_masked

    except Exception as e:
      error =  'IOTA IMPORTER ERROR: Mask failed! {}'.format(e)
      print (error)

    return data, error

  def rename_converted_image(self, filepath):
    # Generate converted image pickle filename
    rename_choice = str(
      self.params.cctbx_ha14.image_conversion.rename_pickle).lower()
    if rename_choice in ("keep_file_structure", "none"):
      img_dir = util.make_image_path(filepath, self.input_base, self.conv_base)
      img_filename = util.make_filename(filepath, new_ext='pickle')
    else:
      self.input_base = self.conv_base
      if rename_choice == "auto_filename":
        prefix = self.init.user_id
      elif rename_choice == "custom_filename":
        prefix = self.params.cctbx_ha14.image_conversion.rename_pickle_prefix
      else:
        prefix = 'iota'
      img_dir = self.conv_base
      number = int(os.path.basename(self.conv_base))
      img_filename = "{}_{}_{:05d}.pickle" \
                 "".format(prefix, number, self.img_object.img_index)
    return os.path.abspath(os.path.join(img_dir, img_filename))

  def apply_threshold(self, data):
    """ Identifies beamstop shadow and sets pixels inside to -2 (to be ignored by
        processing software). Clunky now (merely finds all pixels with intensity
        less than 0.4 * image average intensity), to be refined later.
    """
    try:
      img_raw_bytes = data['DATA']
      beamstop = self.params.image_import.beamstop

      img_thresh = int((beamstop*flex.mean(img_raw_bytes.as_double())))
      beam_stop_sel = img_raw_bytes <= img_thresh
      img_masked = img_raw_bytes.set_selected(beam_stop_sel, -2)

      # mask extensive overloads, too
      #top_thresh = data['SATURATED_VALUE']
      #beam_stop_sel = img_raw_bytes >= top_thresh
      #img_masked_2 = img_masked.set_selected(beam_stop_sel, -2)

      data['DATA'] = img_masked
      return data, None

    except Exception as e:
      error = 'IOTA IMPORTER ERROR: Image thresholding failed! {}'.format(e)
      return None, error


  def square_pickle(self, data):
    """ A function to crop the image pickle to a square such that the beam center
        is in the center of image (mostly copied from cxi_image2pickle.py)
        CCTBX.XFEL ONLY """

    error = None

    # only one active area is allowed, and it should be the size of the image.
    try:
      test = flex.int([0, 0, data['SIZE1'], data['SIZE2']]) == data['ACTIVE_AREAS']
      assert test.all_eq(True)
    except AssertionError as e:
      error = 'IOTA IMPORTER ERROR: Image modification failed! {}'.format(e)
      return None, error

    # retrieve parameters from the dictionary
    pixel_size = data['PIXEL_SIZE']
    beam_x = int(round(data['BEAM_CENTER_X'] / pixel_size))
    beam_y = int(round(data['BEAM_CENTER_Y'] / pixel_size))
    width = data['SIZE1']
    height = data['SIZE2']
    pixels = data['DATA']
    left = beam_x
    right = width - left
    top = beam_y
    bottom = height - top
    new_pixels = pixels.deep_copy()
    new_size = None

    # the new image will be twice the size as the smallest distance between the
    # beam center and one of the edges of the images
    if self.params.cctbx_ha14.image_conversion.square_mode == 'crop':
      try:
        new_half_size = min([right, left, top, bottom])
        new_size = new_half_size * 2
        min_x = int(beam_x - new_half_size)
        min_y = int(beam_y - new_half_size)
        max_x = int(beam_x + new_half_size)
        max_y = int(beam_y + new_half_size)

        if self.params.advanced.flip_beamXY:
          new_beam_x = data['BEAM_CENTER_X'] - min_y * pixel_size
          new_beam_y = data['BEAM_CENTER_Y'] - min_x * pixel_size
          new_pixels = pixels[min_x:max_x, min_y:max_y]
        else:
          new_beam_x = data['BEAM_CENTER_X'] - min_x * pixel_size
          new_beam_y = data['BEAM_CENTER_Y'] - min_y * pixel_size
          new_pixels = pixels[min_y:max_y, min_x:max_x]

        assert new_pixels.focus()[0] == new_pixels.focus()[1]
      except (AssertionError, Exception) as e:
        error = 'IOTA IMPORTER ERROR: Image modification failed! {}'.format(e)
        return None, error

    # the new image will be twice the size as the LARGEST distance between the
    # beam center and one of the edges of the images
    elif self.params.cctbx_ha14.image_conversion.square_mode == 'pad':
      new_half_size = max([right, left, top, bottom])
      new_size = new_half_size * 2
      delta_x = new_half_size - beam_x
      delta_y = new_half_size - beam_y
      new_pixels.resize(flex.grid(new_size, new_size))
      new_pixels.fill(-2)

      if self.params.advanced.flip_beamXY:
        new_beam_x = data['BEAM_CENTER_X'] + delta_y * pixel_size
        new_beam_y = data['BEAM_CENTER_Y'] + delta_x * pixel_size
        new_pixels.matrix_paste_block_in_place(pixels, delta_x, delta_y)
      else:
        new_beam_x = data['BEAM_CENTER_X'] + delta_x * pixel_size
        new_beam_y = data['BEAM_CENTER_Y'] + delta_y * pixel_size
        new_pixels.matrix_paste_block_in_place(pixels, delta_y, delta_x)

    else:
      new_beam_x = data['BEAM_CENTER_X']
      new_beam_y = data['BEAM_CENTER_Y']

    # save the results
    if new_size is not None:
      data['DATA'] = new_pixels
      data['SIZE1'] = new_size
      data['SIZE2'] = new_size
      data['ACTIVE_AREAS'] = flex.int([0, 0, new_size, new_size])
      data['BEAM_CENTER_X'] = new_beam_x
      data['BEAM_CENTER_Y'] = new_beam_y

    return data, error

  def modify_image(self, data=None):
    if not data:
      error = 'IOTA IMPORTER ERROR: Modification failed, no data!'
      return None, error

    # Deactivate image squaring if beam center is in the center of the image
    # CCTBX.XFEL ONLY (Need a better displacement cutoff)
    if (abs(data['BEAM_CENTER_X'] - data['BEAM_CENTER_Y']) < 0.1) :
      self.params.cctbx_ha14.image_conversion.square_mode = 'no_modification'


    # Check if conversion/modification is required and carry them out
    if (
            self.params.cctbx_ha14.image_conversion.square_mode != "no_modification" or
            self.params.image_import.beam_center.x != 0 or
            self.params.image_import.beam_center.y != 0 or
            self.params.image_import.beamstop != 0 or
            self.params.image_import.distance != 0
    ):

      # Convert raw image to image pickle
      beamstop = self.params.image_import.beamstop
      distance = self.params.image_import.distance
      beam_center = [self.params.image_import.beam_center.x,
                     self.params.image_import.beam_center.y]
      square = self.params.cctbx_ha14.image_conversion.square_mode
      mask_file = self.params.image_import.mask

      # Apply mask
      if mask_file is not None:
        data, error = self.apply_mask_from_file(data, mask_file)
        if not data:
          return None, error

      # Override beam center
      if beam_center != [0,0]:
        pixel_size = data['PIXEL_SIZE']
        data['BEAM_CENTER_X'] = int(round(beam_center[0] * pixel_size))
        data['BEAM_CENTER_Y'] = int(round(beam_center[1] * pixel_size))
      if distance != 0:
        data['DISTANCE'] = distance
      if square in ('crop', 'pad'):
        data, error = self.square_pickle(data)
        if not data:
          return None, error
      if beamstop != 0:
        data, error = self.apply_threshold(data)
        if not data:
          return None, error

    # Generate new name for a converted image pickle
    new_path = self.rename_converted_image(self.img_object.img_path)
    self.img_object.raw_img = self.img_object.img_path
    self.img_object.conv_img = new_path
    self.img_object.img_path = new_path

    # Output converted image to file
    try:
      img_dir = os.path.dirname(new_path)
      if not os.path.isdir(img_dir):
        os.makedirs(img_dir)
    except OSError as e:
      error = 'IOTA IMPORTER ERROR: Cannot create folder! {}'.format(e)
      return error, None
    else:
      ep.dump(new_path, data)

    return data, None

  def make_image_object(self, input_entry):
    self.img_object, error = self.import_image(input_entry)

    # Triage image and generate grid search parameters
    if self.img_object.datablock:
      triage_type = str(self.params.cctbx_ha14.image_triage.type).lower()
      if triage_type in ('simple', 'grid_search'):
        hmed, amed = self.image_triage()
      else:
        hmed = self.params.cctbx_ha14.grid_search.height_median
        amed = self.params.cctbx_ha14.grid_search.area_median
      hr = self.params.cctbx_ha14.grid_search.height_range
      ar = self.params.cctbx_ha14.grid_search.area_range
      sig_search = self.params.cctbx_ha14.grid_search.sig_height_search
      gs_type = str(self.params.cctbx_ha14.grid_search.type).lower()

      # Generate grid in image object (put that in there so that the HA14
      # integrator can work)
      if not self.img_object.fail:
        self.img_object.generate_grid(gs_type=gs_type, sig_search=sig_search,
                                      hr=hr, ar=ar, hmed=hmed, amed=amed)
        self.img_object.status = 'imported'
      else:
        self.img_object.status = 'final'

    return self.img_object

class SingleImage(SingleImageBase):    # For current cctbx.xfel
  def __init__(self, imgpath, idx=None):
    SingleImageBase.__init__(self, imgpath=imgpath, idx=idx)
    self.center_int = None
    self.gain = 1.0
    self.img_index = idx

class ImageImporter(ImageImporterBase):
  def __init__(self, init):
    ImageImporterBase.__init__(self, init=init)

  def instantiate_image_object(self, filepath, idx=None):
    ''' Instantiate a SingleImage object for current backend
    :param filepath: path to image file
    '''
    self.img_object = SingleImage(imgpath=filepath, idx=idx)

  def calculate_parameters(self, datablock=None):
    ''' Image modification for current cctbx.xfel '''

    if not datablock:
      # If data are given, apply modifications as specified below
      error = 'IOTA IMPORT ERROR: Datablock not found!'
      return None, error
    else:
      error = []

      # Calculate auto-threshold
      # TODO: Revisit this; I REALLY don't like it.
      if self.iparams.cctbx_xfel.auto_threshold:
        beamX = self.img_object.final['beamX']
        beamY = self.img_object.final['beamY']
        px_size = self.img_object.final['pixel']
        try:
          img = dxtbx.load(self.img_object.img_path)
          raw_data = img.get_raw_data()
          beam_x_px = int(beamX / px_size)
          beam_y_px = int(beamY / px_size)
          data_array = raw_data.as_numpy_array().astype(float)
          self.center_int = np.nanmax(data_array[beam_y_px - 20:beam_y_px + 20,
                                                 beam_x_px - 20:beam_x_px + 20])
        except Exception as e:
          error.append('IOTA IMPORT ERROR: Auto-threshold failed! {}'.format(e))

      # Estimate gain (or set gain to 1.00 if cannot calculate)
      if self.iparams.advanced.estimate_gain:
        from dials.command_line.estimate_gain import estimate_gain
        with util.Capturing() as junk_output:
          try:
            assert self.img_object.datablock   # Must have datablock here
            imageset = self.img_object.datablock.extract_imagesets()[0]
            self.img_object.gain = estimate_gain(imageset)
          except Exception as e:
            error.append('IOTA IMPORT ERROR: Estimate gain failed! '.format(e))

      # Collect error messages for logging
      if error:
        error_message = '\n'.join(error)
      else:
        error_message = None

      return datablock, error_message


# **************************************************************************** #
