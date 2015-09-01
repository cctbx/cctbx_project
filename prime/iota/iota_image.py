from __future__ import division

'''
Author      : Lyubimov, A.Y.
Created     : 10/10/2014
Last Changed: 08/31/2015
Description : Creates image object. If necessary, converts raw image to pickle
              files; crops or pads pickle to place beam center into center of
              image; masks out beam stop. (Adapted in part from
              cxi_image2pickle.py by Aaron Brewster.) If selected, checks image
              for diffraction. Creates instances of integrator and selector
              objects.
'''

import os
import math

from scitbx.array_family import flex

import dxtbx
from cPickle import load
from libtbx import easy_pickle as ep
from xfel.cxi.cspad_ana.cspad_tbx import dpack, evt_timestamp

import prime.iota.iota_misc as misc
import prime.iota.iota_vis_integration as viz

class Empty: pass

class SingleImage(object):

  def __init__(self, img, init, verbose=True, imported_grid=None):
    """ Constructor for the SingleImage object using a raw image file or pickle
    """

    # Initialize parameters
    self.params = init.params
    self.args = init.args
    self.raw_img = img[2]
    self.conv_img = img[2]
    self.img_index = img[0]
    self.img_type = None
    self.status = None
    self.fail = None
    self.Bragg = 0
    self.log_info = []
    self.gs_results = []
    self.main_log = init.logfile
    self.verbose = verbose

    self.input_base = init.input_base
    self.conv_base = init.conv_base
    self.int_base = init.int_base
    self.gs_base = init.gs_base
    self.fin_base = init.fin_base
    self.viz_base = init.viz_base

    self.gs_path = None
    self.gs_file = None
    self.fin_path = None
    self.fin_file = None
    self.viz_path = None

    if not self.params.grid_search.flag_on:
      self.params.grid_search.area_range = 0
      self.params.grid_search.height_range = 0
      self.params.grid_search.sig_height_search = False
    self.grid, self.final = self.generate_grid()

  def import_image(self):
    """ Returns imported image object """
    self.status = 'imported'
    return self


  def import_int_file(self, init):
    """ Replaces path settings in imported image object with new settings
        NEED TO RE-DO LATER """

    # Generate paths to output files
    self.params = init.params
    self.main_log = init.logfile
    self.input_base = init.input_base
    self.conv_base = init.conv_base
    self.int_base = init.int_base
    self.gs_base = init.gs_base
    self.fin_base = init.fin_base
    self.viz_base = init.viz_base
    self.gs_path = misc.make_image_path(self.conv_img, self.input_base, self.gs_base)
    self.gs_file = os.path.abspath(os.path.join(self.gs_path,
            os.path.basename(self.conv_img).split('.')[0] + ".int"))
    self.fin_path = misc.make_image_path(self.conv_img, self.input_base, self.fin_base)
    self.fin_file = os.path.abspath(os.path.join(self.fin_path,
            os.path.basename(self.conv_img).split('.')[0] + "_int.pickle"))
    self.final['final'] = self.fin_file
    self.final['img'] = self.conv_img
    if self.viz_base != None:
      self.viz_path = misc.make_image_path(self.raw_img, self.input_base, self.viz_base)
      self.viz_file = os.path.join(self.viz_path,
                   os.path.basename(self.conv_img).split('.')[0] + "_int.png")

    # Generate output folders and files:
    # Grid search subfolder or final integration subfolder
    if not os.path.isdir(self.gs_path):
      os.makedirs(self.gs_path)
    if not os.path.isdir(self.fin_path):
      os.makedirs(self.fin_path)

    # Grid search / integration log file
    self.int_log = os.path.join(self.fin_path,
                          os.path.basename(self.conv_img).split('.')[0] + '.log')

    # Visualization subfolder
    if self.viz_path != None and not os.path.isdir(self.viz_path):
        os.makedirs(self.viz_path)

    # Reset status to 'grid search' to pick up at selection (if no fail)
    if self.fail == None:
      self.status = 'grid search'

    return self


  def generate_grid(self):
    """ Function to generate grid search parameters for this image object """

    gs_block = []
    h_min = self.params.grid_search.height_median - self.params.grid_search.height_range
    h_max = self.params.grid_search.height_median + self.params.grid_search.height_range
    h_std = self.params.grid_search.height_range
    a_min = self.params.grid_search.area_median - self.params.grid_search.area_range
    a_max = self.params.grid_search.area_median + self.params.grid_search.area_range
    a_std = self.params.grid_search.area_range

    for spot_area in range(a_min, a_max + 1):
      for spot_height in range (h_min, h_max + 1):
        if self.params.grid_search.sig_height_search:
          if spot_height >= 1 + h_std:
            sigs = range(spot_height - h_std, spot_height + 1)
          elif spot_height < 1 + h_std:
            sigs = range(1, spot_height + 1)
          elif spot_height == 1:
            sigs = [1]
        else:
          sigs = [spot_height]

        for sig_height in sigs:
          gs_item = {'sih':sig_height, 'sph':spot_height,
                     'spa':spot_area, 'a':0, 'b':0, 'c':0,
                     'alpha':0, 'beta':0, 'gamma':0, 'sg':'',
                     'strong':0, 'res':0, 'mos':0,
                     'epv':0, 'info':'', 'ok':True}
          gs_block.append(gs_item)

    int_line = {'img':self.conv_img, 'sih':0, 'sph':0, 'spa':0, 'a':0, 'b':0,
                'c':0, 'alpha':0, 'beta':0, 'gamma':0, 'sg':'', 'strong':0,
                'res':0, 'mos':0, 'epv':0, 'info':'', 'final':None}

    return gs_block, int_line


  def load_image(self):
    """ Reads raw image file and extracts data for conversion into pickle
        format. """
    try:
      with misc.Capturing() as junk_output:
        loaded_img = dxtbx.load(self.raw_img)
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
        if abs(beam_x - beam_y) <= 0.1 or self.params.image_conversion.square_mode == "None":
          img_type = 'converted'
        else:
          img_type = 'unconverted'
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
        img_type = 'unconverted'
        if osc_start != osc_range:
          data['OSC_START'] = osc_start
          data['OSC_RANGE'] = osc_range
          data['TIME'] = scan.get_exposure_times()[0]
    else:
      data = None
    return data, img_type


  def square_pickle(self, data):
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
    if self.params.image_conversion.square_mode == 'crop':
      new_half_size = min([right, left, top, bottom])
      new_size = new_half_size * 2
      min_x = beam_x - new_half_size
      min_y = beam_y - new_half_size
      new_beam_x = data['BEAM_CENTER_X'] - min_x * pixel_size
      new_beam_y = data['BEAM_CENTER_Y'] - min_y * pixel_size
      new_pixels = pixels[min_y:min_y+new_size,min_x:min_x+new_size]
      assert new_pixels.focus()[0] == new_pixels.focus()[1]

    # the new image will be twice the size as the LARGEST distance between the
    # beam center and one of the edges of the images
    elif self.params.image_conversion.square_mode == 'pad':
      new_half_size = max([right, left, top, bottom])
      new_size = new_half_size * 2
      min_x = 0
      min_y = 0
      delta_x = new_half_size - beam_x
      delta_y = new_half_size - beam_y
      new_beam_x = data['BEAM_CENTER_X'] + delta_x * pixel_size
      new_beam_y = data['BEAM_CENTER_Y'] + delta_y * pixel_size

      new_pixels.resize(flex.grid(new_size, new_size))
      new_pixels.fill(-2)

      new_pixels.matrix_paste_block_in_place(pixels, delta_y, delta_x)


    # save the results
    data['DATA'] = new_pixels
    data['SIZE1'] = new_size
    data['SIZE2'] = new_size
    data['BEAM_CENTER_X'] = new_beam_x
    data['BEAM_CENTER_Y'] = new_beam_y
    data['ACTIVE_AREAS'] = flex.int([0 , 0, new_size, new_size])

    return data


  def mask_image(self, data):
    """ Identifies beamstop shadow and sets pixels inside to -2 (to be ignored by
        processing software). Clunky now (merely finds all pixels with intensity
        less than 0.4 * image average intensity), to be refined later.

        input: data - image data

        output: data - modified image data
    """
    img_raw_bytes = data['DATA']
    beamstop = self.params.image_conversion.beamstop

    img_thresh = int((beamstop*flex.mean(img_raw_bytes.as_double())))
    beam_stop_sel = img_raw_bytes <= img_thresh
    img_masked = img_raw_bytes.set_selected(beam_stop_sel, -2)

    data['DATA'] = img_masked
    return data


  def convert_image(self):
    """ Converts images into pickle format; crops and masks out beamstop if
        selected """

    img_data, self.img_type = self.load_image()
    info = []

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


    if self.img_type == 'unconverted':
      # Check for and/or create a converted pickles folder
      self.input_base = self.conv_base
      if not os.path.isdir(self.conv_base):
        os.makedirs(self.conv_base)

      # Generate converted image pickle filename
      if self.params.image_conversion.rename_pickle_prefix != None:
        if str(self.params.image_conversion.rename_pickle_prefix).lower() == "auto":
          prefix = os.getlogin()
          suffix = int(os.path.basename(self.conv_base))
        else:
          prefix = self.params.image_conversion.rename_pickle_prefix
          suffix = int(os.path.basename(self.conv_base))
        self.conv_img = os.path.abspath(os.path.join(self.conv_base,
                    "{}_{}_{:05d}.pickle".format(prefix, suffix, self.img_index)))
      else:
        img_path = misc.make_image_path(self.raw_img, self.input_base, self.conv_base)
        self.conv_img = os.path.abspath(os.path.join(img_path,
             os.path.basename(self.raw_img).split('.')[0] + "_converted.pickle"))

      # Convert raw image to image pickle
      beamstop = self.params.image_conversion.beamstop
      beam_center = [self.params.image_conversion.beam_center.x,
                     self.params.image_conversion.beam_center.y]
      square = self.params.image_conversion.square_mode
      if beam_center != [0,0]:
        pixel_size = img_data['PIXEL_SIZE']
        img_data['BEAM_CENTER_X'] = int(round(beam_center[0] * pixel_size))
        img_data['BEAM_CENTER_Y'] = int(round(beam_center[1] * pixel_size))
      if square != "None":
        img_data = self.square_pickle(img_data)
      if beamstop != 0:
        img_data = self.mask_image(img_data)
      self.log_info.append('Converted image : {}'.format(self.conv_img))
      self.log_info.append('Parameters      : BEAM_X = {:<4.2f}, BEAM_Y = {:<4.2f}, '\
                           'PIXEL_SIZE = {:<8.6f}, IMG_SIZE = {:<4} X {:<4}, '\
                           'DIST = {}'.format(img_data['BEAM_CENTER_X'],
                                              img_data['BEAM_CENTER_Y'],
                                              img_data['PIXEL_SIZE'],
                                              img_data['SIZE1'],
                                              img_data['SIZE2'],
                                              img_data['DISTANCE']))
      self.img_type = 'converted'


      # Save converted image pickle
      ep.dump(self.conv_img, img_data)

    self.status = 'converted'
    if self.params.image_triage.flag_on and self.status == 'converted':
      self.fail = self.triage_image()
      self.status = 'triaged'
    else:
      self.fail = None

    if not self.params.image_conversion.convert_only:
      self.gs_path = misc.make_image_path(self.conv_img, self.input_base, self.gs_base)
      self.gs_file = os.path.abspath(os.path.join(self.gs_path,
              os.path.basename(self.conv_img).split('.')[0] + ".int"))
      self.fin_path = misc.make_image_path(self.conv_img, self.input_base, self.fin_base)
      self.fin_file = os.path.abspath(os.path.join(self.fin_path,
              os.path.basename(self.conv_img).split('.')[0] + "_int.pickle"))
      self.final['final'] = self.fin_file
      self.final['img'] = self.conv_img
      if self.viz_base != None:
        self.viz_path = misc.make_image_path(self.raw_img, self.input_base, self.viz_base)
        self.viz_file = os.path.join(self.viz_path,
                   os.path.basename(self.conv_img).split('.')[0] + "_int.png")


      # Generate output folders and files:
      # Grid search subfolder or final integration subfolder
      if not os.path.isdir(self.gs_path):
        os.makedirs(self.gs_path)
      if not os.path.isdir(self.fin_path):
        os.makedirs(self.fin_path)

      # Grid search / integration log file
      self.int_log = os.path.join(self.fin_path,
                            os.path.basename(self.conv_img).split('.')[0] + '.log')

      # Visualization subfolder
      if self.viz_path != None and not os.path.isdir(self.viz_path):
          os.makedirs(self.viz_path)

    return self


  def triage_image(self):
    """ Performs a quick DISTL spotfinding without grid search.
        NOTE: Convert to DIALS spotfinder when ready! (Even for CCTBX runs.)
    """

    import spotfinder
    from spotfinder.command_line.signal_strength import master_params as sf_params
    from spotfinder.applications.wrappers import DistlOrganizer

    sf_params = sf_params.extract()
    sf_params.distl.image = self.conv_img

    E = Empty()
    E.argv=['Empty']
    E.argv.append(sf_params.distl.image)

    selected_output = []
    total_output = []
    bragg_spots = []
    spotfinding_log = ['{}\n'.format(self.conv_img)]

    # set spotfinding parameters for DISTL spotfinder
    sf_params.distl.minimum_spot_area = self.params.grid_search.area_median
    sf_params.distl.minimum_spot_height = self.params.grid_search.height_median
    sf_params.distl.minimum_signal_height = int(self.params.grid_search.height_median / 2)

    # run DISTL spotfinder
    with misc.Capturing() as junk_output:
      Org = DistlOrganizer(verbose = False, argument_module=E,
                           phil_params=sf_params)

      Org.printSpots()

    # Extract relevant spotfinding info & make selection
    for frame in Org.S.images.keys():
      self.Bragg = Org.S.images[frame]['N_spots_inlier']

    if self.Bragg >= self.params.image_triage.min_Bragg_peaks:
      self.log_info.append('ACCEPTED! {} good Bragg peaks'\
                            ''.format(self.Bragg))
      status = None
    else:
      self.log_info.append('REJECTED')
      status = 'failed triage'

    return status


  def determine_gs_result_file(self):
    if self.params.selection.select_only.grid_search_path != None:
      gs_path = os.path.abspath(self.params.selection.select_only.grid_search_path)
    else:
      run_number = int(os.path.basename(self.int_base)) - 1
      gs_path = "{}/integration/{:03d}/grid_search"\
                "".format(os.path.abspath(os.curdir), run_number)
    gs_result_file = os.path.join(gs_path, os.path.basename(self.gs_file))
    return gs_result_file


  def integrate_cctbx(self, tag):
    """ Runs integration using the Integrator class """

    # Check to see if the image is suitable for grid search / integration
    if self.fail != None:
      self.grid = []
      self.final['final'] = None
    else:
      from prime.iota.iota_integrate import Integrator
      integrator = Integrator(self.conv_img,
                              self.fin_file,
                              self.params.selection.min_sigma,
                              self.params.target,
                              self.params.advanced.charts,
                              self.viz_path,
                              self.int_log,
                              tag)
      if tag == 'grid search':
        self.log_info.append('\nCCTBX grid search:')
        for i in range(len(self.grid)):
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

      elif tag == 'integrate':
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

        if self.params.advanced.viz == 'integration':
          viz.make_png(self.final['img'], self.final['final'], self.viz_file)
        elif self.params.advanced.viz == 'cv_vectors':
          viz.cv_png(self.final['img'], self.final['final'], self.viz_file)

      # Save image object to file
      ep.dump(self.gs_file, self)

    return self


  def select_cctbx(self):
    """ Selects best grid search result using the Selector class """

    if self.fail == None:
      from prime.iota.iota_integrate import Selector
      selector = Selector(self.grid,
                          self.final,
                          self.params.selection.prefilter.flag_on,
                          self.params.selection.prefilter.target_uc_tolerance,
                          self.params.selection.prefilter.target_pointgroup,
                          self.params.selection.prefilter.target_unit_cell,
                          self.params.selection.prefilter.min_reflections,
                          self.params.selection.prefilter.min_resolution)

      self.fail, self.final, log_entry = selector.select()
      self.status = 'selection'
      self.log_info.append(log_entry)

    # Save results into a pickle file
    ep.dump(self.gs_file, self)

    return self

  def process(self):
    """ Image processing; selects method, runs requisite modules """

    if self.params.advanced.integrate_with == 'cctbx':

      # Run grid search if haven't already
      if self.fail == None and self.status != 'grid search':
        self.integrate_cctbx('grid search')

      # Run selection if haven't already
      if self.fail == None and self.status != 'select':
        self.select_cctbx()

      # Run final integration if haven't already
      if self.fail == None and self.status != 'final':
        self.integrate_cctbx('integrate')

      if self.verbose:
         log_entry = "\n".join(self.log_info)
         misc.main_log(self.main_log, log_entry)
         misc.main_log(self.main_log, '\n{:-^100}\n'.format(''))

    elif self.params.advanced.integrate_with == 'dials':
      print "PLACEHOLDER FOR DIALS INTEGRATION!"

    return self


# **************************************************************************** #
