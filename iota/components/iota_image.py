from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 02/12/2018
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

from scitbx.array_family import flex

import dxtbx
from libtbx import easy_pickle as ep
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp

import iota.components.iota_misc as misc
import iota.components.iota_vis_integration as viz

class SingleImage(object):

  def __init__(self, img, init, verbose=True, imported_grid=None):
    """ Constructor for the SingleImage object using a raw image file or pickle
    """

    # Initialize parameters
    self.params = init.params
    self.args = init.args
    self.user_id = init.user_id
    self.raw_img = img[2]
    self.conv_img = img[2]
    self.img_index = img[0]
    self.center_int = 0
    self.status = None
    self.fail = None
    self.final = None
    self.log_info = []
    self.gs_results = []
    self.main_log = init.logfile
    self.verbose = verbose
    self.hmed = self.params.cctbx.grid_search.height_median
    self.amed = self.params.cctbx.grid_search.area_median

    self.input_base = init.input_base
    self.conv_base = init.conv_base
    self.int_base = init.int_base
    self.obj_base = init.obj_base
    self.fin_base = init.fin_base
    self.viz_base = init.viz_base
    self.log_base = init.log_base
    self.tmp_base = init.tmp_base
    self.abort_file = os.path.join(self.int_base, '.abort.tmp')

    self.obj_path = None
    self.obj_file = None
    self.fin_path = None
    self.fin_file = None
    self.viz_path = None


# ============================== SELECTION-ONLY FUNCTIONS ============================== #

  def import_int_file(self, init):
    """ Replaces path settings in imported image object with new settings
        NEED TO RE-DO LATER """

    if os.path.isfile(self.abort_file):
      self.fail = 'aborted'
      return self

    # Generate paths to output files
    self.params = init.params
    self.main_log = init.logfile
    self.input_base = init.input_base
    self.conv_base = init.conv_base
    self.int_base = init.int_base
    self.obj_base = init.obj_base
    self.fin_base = init.fin_base
    self.log_base = init.log_base
    self.viz_base = init.viz_base
    self.obj_path = misc.make_image_path(self.conv_img, self.input_base, self.obj_base)
    self.obj_file = os.path.abspath(os.path.join(self.obj_path,
            os.path.basename(self.conv_img).split('.')[0] + ".int"))
    self.fin_path = misc.make_image_path(self.conv_img, self.input_base, self.fin_base)
    self.log_path = misc.make_image_path(self.conv_img, self.input_base, self.log_base)
    self.fin_file = os.path.abspath(os.path.join(self.fin_path,
            os.path.basename(self.conv_img).split('.')[0] + "_int.pickle"))
    self.final['final'] = self.fin_file
    self.final['img'] = self.conv_img
    self.viz_path = misc.make_image_path(self.conv_img, self.input_base, self.viz_base)
    self.viz_file = os.path.join(self.viz_path,
                    os.path.basename(self.conv_img).split('.')[0] + "_int.png")

    # Create actual folders (if necessary)
    try:
      if not os.path.isdir(self.obj_path):
        os.makedirs(self.obj_path)
      if not os.path.isdir(self.fin_path):
        os.makedirs(self.fin_path)
      if not os.path.isdir(self.viz_path):
        os.makedirs(self.viz_path)
    except OSError:
      pass

    # Grid search / integration log file
    self.int_log = os.path.join(self.log_path,
                          os.path.basename(self.conv_img).split('.')[0] + '.tmp')

    # Reset status to 'grid search' to pick up at selection (if no fail)
    if self.fail == None:
      self.status = 'bypass grid search'

    return self

  def determine_gs_result_file(self):
    """ For 'selection-only' cctbx.xfel runs, determine where the image objects are """
    if self.params.cctbx.selection.select_only.grid_search_path != None:
      obj_path = os.path.abspath(self.params.cctbx.selection.select_only.grid_search_path)
    else:
      run_number = int(os.path.basename(self.int_base)) - 1
      obj_path = "{}/integration/{:03d}/image_objects"\
                "".format(os.path.abspath(os.curdir), run_number)
    gs_result_file = os.path.join(obj_path, os.path.basename(self.obj_file))
    return gs_result_file


# =============================== IMAGE IMPORT FUNCTIONS =============================== #

  def load_image(self):
    """ Reads raw image file and extracts data for conversion into pickle
        format. Also estimates gain if turned on."""
    # Load raw image or image pickle

    try:
      with misc.Capturing() as junk_output:
        loaded_img = dxtbx.load(self.raw_img)
    except Exception, e:
      print 'IOTA IMPORT ERROR:', e
      loaded_img = None
      pass

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
        img_type = 'pickle'
      else:
        img_type = 'raw'
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
    else:
      data = None
      img_type = 'not imported'

    # Estimate gain (or set gain to 1.00 if cannot calculate)
    # Cribbed from estimate_gain.py by Richard Gildea
    if self.params.advanced.estimate_gain:
      from dxtbx.datablock import DataBlockFactory
      from dials.command_line.estimate_gain import estimate_gain
      with misc.Capturing() as junk_output:
        try:
          datablock = DataBlockFactory.from_filenames([self.raw_img])[0]
          imageset = datablock.extract_imagesets()[0]
          self.gain = estimate_gain(imageset)
        except Exception, e:
          self.gain = 1.0
    else:
      self.gain = 1.0

    return data, img_type


  def generate_grid(self):
    """ Function to generate grid search parameters for this image object
        CCTBX.XFEL ONLY """

    self.grid = []
    h_min = self.hmed - self.hrange
    h_max = self.hmed + self.hrange
    a_min = self.amed - self.arange
    a_max = self.amed + self.arange
    h_std = self.params.cctbx.grid_search.height_range
    a_std = self.params.cctbx.grid_search.area_range

    for spot_area in range(a_min, a_max + 1):
      for spot_height in range (h_min, h_max + 1):
        if self.params.cctbx.grid_search.sig_height_search:
          if spot_height >= 1 + h_std:
            sigs = range(spot_height - h_std, spot_height + 1)
          elif spot_height < 1 + h_std:
            sigs = range(1, spot_height + 1)
          elif spot_height == 1:
            sigs = [1]
        else:
          sigs = [spot_height]

        for sig_height in sigs:
          if (spot_area, spot_height, sig_height) not in self.grid_points:
            self.grid_points.append((spot_area, spot_height, sig_height))
            self.grid.append({'sih':sig_height, 'sph':spot_height,
                                  'spa':spot_area, 'a':0, 'b':0, 'c':0,
                                  'alpha':0, 'beta':0, 'gamma':0, 'sg':'',
                                  'strong':0, 'res':0, 'lres':0, 'mos':0,
                                  'epv':0, 'info':'', 'ok':True})

    self.final = {'img':self.conv_img, 'sih':0, 'sph':0, 'spa':0, 'a':0, 'b':0,
                  'c':0, 'alpha':0, 'beta':0, 'gamma':0, 'sg':'', 'strong':0,
                  'res':0, 'lres':0, 'mos':0, 'epv':0, 'info':'', 'final':None,
                  'program':'cctbx'}

  def square_pickle(self, data):
    """ A function to crop the image pickle to a square such that the beam center
        is in the center of image (mostly copied from cxi_image2pickle.py)
        CCTBX.XFEL ONLY """

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
    left = beam_x
    right = width - left
    top = beam_y
    bottom = height - top
    new_pixels = pixels.deep_copy()
    new_size = None

    # the new image will be twice the size as the smallest distance between the
    # beam center and one of the edges of the images
    if self.params.image_conversion.square_mode == 'crop':
      new_half_size = min([right, left, top, bottom])
      new_size = new_half_size * 2
      min_x = beam_x - new_half_size
      min_y = beam_y - new_half_size
      max_x = beam_x + new_half_size
      max_y = beam_y + new_half_size

      if self.params.advanced.flip_beamXY:
        new_beam_x = data['BEAM_CENTER_X'] - min_y * pixel_size
        new_beam_y = data['BEAM_CENTER_Y'] - min_x * pixel_size
        new_pixels = pixels[min_x:max_x, min_y:max_y]
      else:
        new_beam_x = data['BEAM_CENTER_X'] - min_x * pixel_size
        new_beam_y = data['BEAM_CENTER_Y'] - min_y * pixel_size
        new_pixels = pixels[min_y:max_y, min_x:max_x]


      assert new_pixels.focus()[0] == new_pixels.focus()[1]

    # the new image will be twice the size as the LARGEST distance between the
    # beam center and one of the edges of the images
    elif self.params.image_conversion.square_mode == 'pad':
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
      data['ACTIVE_AREAS'] = flex.int([0 , 0, new_size, new_size])
      data['BEAM_CENTER_X'] = new_beam_x
      data['BEAM_CENTER_Y'] = new_beam_y

    return data

  def mask_image(self, data):
    """ Identifies beamstop shadow and sets pixels inside to -2 (to be ignored by
        processing software). Clunky now (merely finds all pixels with intensity
        less than 0.4 * image average intensity), to be refined later.
    """
    img_raw_bytes = data['DATA']
    beamstop = self.params.image_conversion.beamstop

    img_thresh = int((beamstop*flex.mean(img_raw_bytes.as_double())))
    beam_stop_sel = img_raw_bytes <= img_thresh
    img_masked = img_raw_bytes.set_selected(beam_stop_sel, -2)

    # mask extensive overloads, too
    #top_thresh = data['SATURATED_VALUE']
    #beam_stop_sel = img_raw_bytes >= top_thresh
    #img_masked_2 = img_masked.set_selected(beam_stop_sel, -2)

    data['DATA'] = img_masked
    return data

  def apply_mask_from_file(self, data, mask_file, invert=False):
    img_raw_bytes = data['DATA']

    full_mask = ep.load(mask_file)
    assert len(full_mask) == 1
    if invert:
      mask = full_mask[0] == False
    else:
      mask = full_mask[0]

    img_masked = img_raw_bytes.set_selected(mask, -2)

    data['DATA'] = img_masked
    return data

  def import_image(self):
    """ Image conversion:
          - Writes out data in pickle format (cctbx.xfel only)
          - Moves beam center into center of image, crops / pads image (cctbx.xfel only)
          - Adjusts beam center and distance (optional)
          - Thresholds and masks beamstop shadow (optional)
    """

    if os.path.isfile(self.abort_file):
      self.fail = 'aborted'
      return self

    # Load image
    img_data, img_type = self.load_image()
    self.status = 'loaded'

    if img_data is None:
      self.log_path = misc.make_image_path(self.conv_img, self.input_base,
                                           self.log_base)
      self.int_log = os.path.join(self.log_path,
                                  os.path.basename(self.conv_img).split('.')[
                                    0] + '.tmp')
      self.log_info.append('\n{:-^100}\n'.format(self.raw_img))
      self.log_info.append('FAILED TO IMPORT')
      self.status = 'failed import'
      self.fail = 'failed import'

      if not self.params.image_conversion.convert_only:
        self.obj_path = misc.make_image_path(self.conv_img, self.input_base,
                                             self.obj_base)
        rej_name = filter(None, os.path.basename(self.conv_img).split('.'))[
                     0] + '_reject.int'
        self.obj_file = os.path.abspath(os.path.join(self.obj_path, rej_name))

        try:
          if not os.path.isdir(self.obj_path):
            os.makedirs(self.obj_path)
        except OSError:
          pass

        ep.dump(self.obj_file, self)

      return self

    # if DIALS is selected, change image type to skip conversion step
    if self.params.advanced.integrate_with == 'dials':
      img_type = 'dials_input'
      if self.params.dials.auto_threshold:
        beam_x_px = img_data['BEAM_CENTER_X'] / img_data['PIXEL_SIZE']
        beam_y_px = img_data['BEAM_CENTER_Y'] / img_data['PIXEL_SIZE']
        data_array = img_data['DATA'].as_numpy_array().astype(float)
        self.center_int = np.nanmax(data_array[beam_y_px - 20:beam_y_px + 20,
                                               beam_x_px - 20:beam_x_px + 20])

    # Log initial image information
    self.log_info.append('\n{:-^100}\n'.format(self.raw_img))
    self.log_info.append('Imported image  : {}'.format(self.raw_img))
    self.log_info.append('Parameters      : BEAM_X = {:<4.2f}, BEAM_Y = {:<4.2f}, '\
                         'PIXEL_SIZE = {:<8.6f}, IMG_SIZE = {:<4} X {:<4}, '\
                         'DIST = {}'.format(img_data['BEAM_CENTER_X'],
                                            img_data['BEAM_CENTER_Y'],
                                            img_data['PIXEL_SIZE'],
                                            img_data['SIZE1'],
                                            img_data['SIZE2'],
                                            img_data['DISTANCE']))

    # Deactivate image squaring if beam center is in the center of the image
    # CCTBX.XFEL ONLY (Need a better displacement cutoff)
    if (                  self.params.advanced.integrate_with == 'dials' or
        abs(img_data['BEAM_CENTER_X'] - img_data['BEAM_CENTER_Y']) < 0.1
        ):
      self.params.image_conversion.square_mode = 'no_modification'

    # Check if conversion/modification is required and carry them out
    if self.params.advanced.integrate_with != 'dials' and \
       (
                                                       img_type == 'raw' or
           self.params.image_conversion.square_mode != "no_modification" or
                         self.params.image_conversion.beam_center.x != 0 or
                         self.params.image_conversion.beam_center.y != 0 or
                              self.params.image_conversion.beamstop != 0 or
                              self.params.image_conversion.distance != 0
        ):

      # Check for and/or create a converted pickles folder
      try:
        if not os.path.isdir(self.conv_base):
          os.makedirs(self.conv_base)
      except OSError:
        pass

      # Generate converted image pickle filename
      rename_choice = str(self.params.image_conversion.rename_pickle).lower()
      if rename_choice in ("keep_file_structure", "none"):
        img_path = misc.make_image_path(self.raw_img, self.input_base,
                                        self.conv_base)
        img_filename = os.path.basename(self.raw_img).split('.')[0] + ".pickle"
        self.conv_img = os.path.abspath(os.path.join(img_path, img_filename))
        try:
          if not os.path.isdir(img_path):
            os.makedirs(img_path)
        except OSError:
          pass
      else:
        if rename_choice == "auto_filename":
          prefix = self.user_id
        elif rename_choice == "custom_filename":
          prefix = self.params.image_conversion.rename_pickle_prefix
        number = int(os.path.basename(self.conv_base))
        self.conv_img = os.path.abspath(os.path.join(self.conv_base,
                    "{}_{}_{:05d}.pickle".format(prefix, number, self.img_index)))

      # Convert raw image to image pickle
      beamstop = self.params.image_conversion.beamstop
      distance = self.params.image_conversion.distance
      beam_center = [self.params.image_conversion.beam_center.x,
                     self.params.image_conversion.beam_center.y]
      square = self.params.image_conversion.square_mode
      mask_file = self.params.image_conversion.mask
      invert = self.params.image_conversion.invert_boolean_mask
      if mask_file is not None:
        img_data = self.apply_mask_from_file(img_data, mask_file, invert)
      if beam_center != [0,0]:
        pixel_size = img_data['PIXEL_SIZE']
        img_data['BEAM_CENTER_X'] = int(round(beam_center[0] * pixel_size))
        img_data['BEAM_CENTER_Y'] = int(round(beam_center[1] * pixel_size))
      if distance != 0:
        img_data['DISTANCE'] = distance
      if square in ('crop', 'pad'):
        img_data = self.square_pickle(img_data)
      if beamstop != 0:
        img_data = self.mask_image(img_data)

      # Log converted image information
      self.log_info.append('Converted image : {}'.format(self.conv_img))
      self.log_info.append('Parameters      : BEAM_X = {:<4.2f}, BEAM_Y = {:<4.2f}, '\
                           'PIXEL_SIZE = {:<8.6f}, IMG_SIZE = {:<4} X {:<4}, '\
                           'DIST = {}'.format(img_data['BEAM_CENTER_X'],
                                              img_data['BEAM_CENTER_Y'],
                                              img_data['PIXEL_SIZE'],
                                              img_data['SIZE1'],
                                              img_data['SIZE2'],
                                              img_data['DISTANCE']))
      self.input_base = self.conv_base
      self.status = 'converted'

      # Save converted image pickle
      ep.dump(self.conv_img, img_data)

    # Triage image (i.e. check for usable diffraction, using selected method)
    if str(self.params.image_triage.type).lower() in ('simple', 'grid_search'):
      if self.params.advanced.integrate_with == 'cctbx':
        from iota.components.iota_cctbx import Triage
        triage = Triage(self.conv_img, self.gain, self.params)
        self.fail, log_entry, self.hmed, self.amed = triage.triage_image()

      elif self.params.advanced.integrate_with == 'dials':
        from iota.components.iota_dials import Triage
        triage = Triage(img=self.conv_img,
                        gain=self.gain,
                        center_intensity=self.center_int,
                        params=self.params)
        self.fail, log_entry = triage.triage_image()

      self.log_info.append(log_entry)
      self.status = 'triaged'
    else:
      self.fail = None

    # Generate integration result dictionary for cctbx.xfel or DIALS
    if self.params.advanced.integrate_with == 'cctbx':
      gs_type = str(self.params.cctbx.grid_search.type).lower()
      if gs_type == 'smart':
        self.hrange = 1
        self.arange = 1
      elif gs_type in ('none', 'no_grid_search'):
        self.hrange = 0
        self.arange = 0
      else:
        self.hrange = self.params.cctbx.grid_search.height_range
        self.arange = self.params.cctbx.grid_search.area_range
      self.grid_points = []
      self.generate_grid()
    elif self.params.advanced.integrate_with == 'dials':
      self.final = {'img':self.conv_img, 'a':0, 'b':0, 'c':0, 'alpha':0,
                    'beta':0, 'gamma':0, 'sg':'','strong':0, 'res':0,
                    'lres':0, 'mos':0, 'epv':0, 'info':'','final':None,
                    'program':'dials'}

    # Generate names for output folders and files:
    if not self.params.image_conversion.convert_only:
      self.obj_path = misc.make_image_path(self.conv_img, self.input_base, self.obj_base)
      self.obj_file = os.path.abspath(os.path.join(self.obj_path,
              os.path.basename(self.conv_img).split('.')[0] + ".int"))
      self.fin_path = misc.make_image_path(self.conv_img, self.input_base, self.fin_base)
      self.log_path = misc.make_image_path(self.conv_img, self.input_base, self.log_base)
      self.fin_file = os.path.abspath(os.path.join(self.fin_path,
        "int_{}.pickle".format(os.path.basename(self.conv_img).split('.')[0])))
      self.final['final'] = self.fin_file
      self.final['img'] = self.conv_img
      self.int_log = os.path.join(self.log_path,
                        os.path.basename(self.conv_img).split('.')[0] + '.tmp')
      self.viz_path = misc.make_image_path(self.conv_img, self.input_base, self.viz_base)
      self.viz_file = os.path.join(self.viz_path,
            "int_{}.png".format(os.path.basename(self.conv_img).split('.')[0]))

      # Create actual folders (if necessary)f
      try:
        if not os.path.isdir(self.obj_path):
          os.makedirs(self.obj_path)
        if not os.path.isdir(self.fin_path):
          os.makedirs(self.fin_path)
        if not os.path.isdir(self.log_path):
          os.makedirs(self.log_path)
        if not os.path.isdir(self.viz_path):
          os.makedirs(self.viz_path)
      except OSError:
        pass

      # Save image object to file
      ep.dump(self.obj_file, self)

    self.status = 'imported'

    # If conversion only option is selected, write conversion info to log
    if self.params.image_conversion.convert_only:
      log_entry = "\n".join(self.log_info)
      misc.main_log(self.main_log, log_entry)

    return self


# ================================= IMAGE INTEGRATORS ================================== #

  def integrate_cctbx(self, tag, grid_point=0, single_image=False):
    """ Runs integration using the Integrator class """

    # Check to see if the image is suitable for grid search / integration
    if self.fail != None:
      self.grid = []
      self.final['final'] = None
    else:
      from iota.components.iota_cctbx import Integrator
      integrator = Integrator(params = self.params,
                              source_image=self.conv_img,
                              output_image=self.fin_file,
                              viz=self.viz_path,
                              log=self.int_log,
                              tag=tag,
                              tmp_base=self.tmp_base,
                              gain=self.gain,
                              single_image=single_image)
      if tag == 'grid search':
        self.log_info.append('\nCCTBX grid search:')
        for i in range(len(self.grid)):

          # If aborted from GUI
          if os.path.isfile(self.abort_file):
            self.fail = 'aborted'
            return self

          int_results = integrator.integrate(self.grid[i])
          self.grid[i].update(int_results)
          img_filename = os.path.basename(self.conv_img)
          log_entry ='{:<{width}}: S = {:<3} H = {:<3} ' \
                     'A = {:<3} ---> {}'.format(img_filename,
                      self.grid[i]['sih'],
                      self.grid[i]['sph'],
                      self.grid[i]['spa'],
                      self.grid[i]['info'],
                      width = len(img_filename) + 2)
          self.log_info.append(log_entry)
          self.gs_results.append(log_entry)

        # Throw out grid search results that yielded no integration
        self.grid = [i for i in self.grid if "not integrated" not in i['info'] and\
                     "no data recorded" not in i['info']]
        self.status = 'grid search'

      elif tag == 'split grid':

        if os.path.isfile(self.abort_file):
          self.fail = 'aborted'
          return self

        self.log_info.append('\nCCTBX INTEGRATION grid search:')
        int_results = integrator.integrate(self.grid[grid_point])
        self.grid[grid_point].update(int_results)
        img_filename = os.path.basename(self.conv_img)
        log_entry ='{:<{width}}: S = {:<3} H = {:<3} ' \
                   'A = {:<3} ---> {}'.format(img_filename,
                                              self.grid[grid_point]['sih'],
                                              self.grid[grid_point]['sph'],
                                              self.grid[grid_point]['spa'],
                                              self.grid[grid_point]['info'],
                                              width = len(img_filename) + 2)
        self.log_info.append(log_entry)
        self.gs_results.append(log_entry)

      elif tag == 'integrate':

        if os.path.isfile(self.abort_file):
          self.fail = 'aborted'
          return self

        self.log_info.append('\nCCTBX final integration:')
        final_results = integrator.integrate(self.final)
        self.final.update(final_results)
        self.status = 'final'
        img_filename = os.path.basename(self.conv_img)
        log_entry ='{:<{width}}: S = {:<3} H = {:<3} ' \
                   'A = {:<3} ---> {}'.format(img_filename,
                    self.final['sih'],
                    self.final['sph'],
                    self.final['spa'],
                    self.final['info'],
                    width = len(img_filename) + 2)
        self.log_info.append(log_entry)

        if self.params.analysis.viz == 'integration':
          viz.make_png(self.final['img'], self.final['final'], self.viz_file)
        elif self.params.analysis.viz == 'cv_vectors':
          viz.cv_png(self.final['img'], self.final['final'], self.viz_file)

      # Save image object to file
      ep.dump(self.obj_file, self)

    return self


  def select_cctbx(self):
    """ Selects best grid search result using the Selector class """

    if os.path.isfile(self.abort_file):
      self.fail = 'aborted'
      return self

    if self.fail == None:
      from iota.components.iota_cctbx import Selector
      selector = Selector(self.grid,
                          self.final,
                          self.params.cctbx.selection.prefilter.flag_on,
                          self.params.cctbx.selection.prefilter.target_uc_tolerance,
                          self.params.cctbx.selection.prefilter.target_pointgroup,
                          self.params.cctbx.selection.prefilter.target_unit_cell,
                          self.params.cctbx.selection.prefilter.min_reflections,
                          self.params.cctbx.selection.prefilter.min_resolution,
                          self.params.cctbx.selection.select_by)

      self.fail, self.final, log_entry = selector.select()
      self.status = 'selection'
      self.log_info.append(log_entry)

    # Save results into a pickle file
    ep.dump(self.obj_file, self)

    return self

  def process(self, single_image=False):
    """ Image processing; selects method, runs requisite modules """
    if self.status != 'bypass grid search':
      self.status = 'processing'

    #for CCTBX indexing / integration
    if self.params.advanced.integrate_with == 'cctbx':
      terminate = False
      prev_status = self.status
      prev_fail = 'first cycle'
      prev_final = self.final
      prev_epv = 9999

      while not terminate:
        if os.path.isfile(self.abort_file):
          self.fail = 'aborted'
          return self

        # Run grid search if haven't already
        if self.fail == None and 'grid search' not in self.status:
          self.integrate_cctbx('grid search', single_image=single_image)

        # Run selection if haven't already
        if self.fail == None and self.status != 'selection':
          self.select_cctbx()

        # If smart grid search is active run multiple rounds until convergence
        if self.params.cctbx.grid_search.type == 'smart':
          if self.fail == None and self.final['epv'] < prev_epv:
            prev_epv = self.final['epv']
            prev_final = self.final
            prev_status = self.status
            prev_fail = self.fail
            self.hmed = self.final['sph']
            self.amed = self.final['spa']
            self.generate_grid()
            self.final['final'] = self.fin_file
            if len(self.grid) == 0:
              self.final = prev_final
              self.status = prev_status
              self.fail = prev_fail
              terminate = True
              continue
            if self.verbose:
              log_entry = '\nNew starting point: H = {}, A = {}\n'\
                          ''.format(self.hmed, self.amed)
              self.log_info.append(log_entry)
          else:
            if prev_fail != 'first cycle':
              self.final = prev_final
              self.status = prev_status
              self.fail = prev_fail
              if self.verbose:
                log_entry = '\nFinal set of parameters: H = {}, A = {}'\
                            ''.format(self.final['sph'], self.final['spa'])
                self.log_info.append(log_entry)
            terminate = True

        # If brute force grid search is selected run one round
        else:
          terminate = True

      # Run final integration if haven't already
      if self.fail == None and self.status != 'final':
        self.integrate_cctbx('integrate', single_image=single_image)

      # If verbose output selected (default), write to main log
      if self.verbose:
         log_entry = "\n".join(self.log_info)
         misc.main_log(self.main_log, log_entry)
         misc.main_log(self.main_log, '\n\n')

      # Make a temporary process log into a final process log
      if os.path.isfile(self.int_log):
        final_int_log = os.path.join(self.log_path,
                                    os.path.basename(self.int_log).split('.')[0] + '.log')
        os.rename(self.int_log, final_int_log)


    # For DIALS integration (WORK IN PROGRESS)
    elif self.params.advanced.integrate_with == 'dials':

      if os.path.isfile(self.abort_file):
        self.fail = 'aborted'
        return self

      if self.fail is not None:
        self.status = 'final'
        ep.dump(self.obj_file, self)
        return self
      else:
        # Create DIALS integrator object
        from iota.components.iota_dials import Integrator
        integrator = Integrator(source_image=self.conv_img,
                                object_folder=self.obj_path,
                                int_folder=self.int_base,
                                final_filename=self.fin_file,
                                final=self.final,
                                logfile=self.int_log,
                                gain=self.gain,
                                center_intensity=self.center_int,
                                params=self.params)

        # Run DIALS
        self.fail, self.final, int_log = integrator.run()
        if self.fail != None:
          self.status = 'final'
        self.log_info.append(int_log)
        log_entry = "\n".join(self.log_info)
        misc.main_log(self.main_log, log_entry)
        misc.main_log(self.main_log, '\n{:-^100}\n'.format(''))

        # Make a temporary process log into a final process log
        final_int_log = self.int_log.split('.')[0] + ".log"
        os.rename(self.int_log, final_int_log)

    self.status = 'final'
    ep.dump(self.obj_file, self)

# **************************************************************************** #
