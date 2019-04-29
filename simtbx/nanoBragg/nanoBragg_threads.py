from __future__ import division, print_function

'''
Author      : Lyubimov, A.Y.
Created     : 12/21/2017
Last Changed: 12/21/2017
Description : SIMTBX (nanoBragg) GUI Threads and PostEvents
'''

import os
import wx
from threading import Thread

from libtbx import easy_run
from scitbx.array_family import flex

from simtbx.nanoBragg import nanoBragg
from simtbx.nanoBragg import shapetype

from iotbx import pdb
from cctbx.eltbx import henke


# ------------------------------ IMAGE SYNTHESIS ----------------------------- #

tp_EVT_NBDONE = wx.NewEventType()
EVT_NBDONE = wx.PyEventBinder(tp_EVT_NBDONE, 1)

class nanoBraggAllDone(wx.PyCommandEvent):
  ''' Send event when finished all cycles  '''

  def __init__(self, etype, eid, pixels=None):
    wx.PyCommandEvent.__init__(self, etype, eid)
    self.pixels = pixels

  def GetValue(self):
    return self.pixels

class nanoBraggThread(Thread):
  ''' Basic spotfinder (with defaults) that could be used to rapidly analyze
  images as they are collected '''

  def __init__(self, parent, params,
               add_background=True,
               add_noise=True,
               randomize=True,
               pixel_scale=None,
               preview=True):
    Thread.__init__(self)
    self.SIM = nanoBragg(detpixels_slowfast=(1231, 1263),
                         pixel_size_mm=0.344,
                         Ncells_abc=(5, 5, 5))
    self.parent = parent
    self.params = params
    self.add_background = add_background
    self.add_noise = add_noise
    self.randomize = randomize
    self.crystal_size = pixel_scale
    self.preview = preview

  def run(self):
    if self.preview:
      self.generate_image()
    else:
      self.generate_image_set()


  def generate_image_set(self):
    pass

  def generate_image(self):
    img_filename = 'test_image.{}'.format(self.params.image_format)
    self.img_file = os.path.join(self.params.output, img_filename)

    if self.params.reference_coordinates is None:
      self.SIM.default_F = 1000
      print('DEBUG: GENERATING FCALC OF {}'.format(self.SIM.default_F))
    else:
      print('DEBUG: GENERATING FCALC FROM COORDINATES')
      self.SIM.default_F = 0
      sfall = self.fcalc_from_pdb(resolution=1.6, algorithm="fft",
                                  wavelength=self.SIM.wavelength_A)
      self.SIM.Fhkl = sfall

    # fastest option, least realistic
    self.SIM.xtal_shape = shapetype.Tophat

    # oversample = 1 is the fastest setting (speed inversely proportional to
    # square of oversample)
    self.SIM.oversample = 1
    self.SIM.wavelength_A = 1.3
    self.SIM.polarization = 1

    self.SIM.distance_mm = 300

    # default orientation is with a axis down the beam, lets pick a random one
    print('DEBUG: RANDOMIZE ORIENTATION - ', self.randomize)
    if self.randomize:
      self.SIM.randomize_orientation()
    # display randomly-picked missetting angles
    print('SIMTBX: MISSETTING ANGLES: ', self.SIM.missets_deg)
    # or an Arndt-Wonacott A matrix (U*B), same as used by mosflm
    print('SIMTBX: ARNDT-WONACOTT A MATRIX', self.SIM.Amatrix)

    # show all parameters
    print('SIMTBX: SHOWING PARAMETERS:')
    self.SIM.show_params()
    print('\n *** \n')

    # print 'SIMTBX: SETTING BEAM PARAMETERS'
    self.SIM.flux = 1e12
    # assumes round beam
    self.SIM.beamsize_mm = 0.1
    self.SIM.exposure_s = 1
    print('SIMTBX: FLUX = ', self.SIM.flux)
    print('SIMTBX: BEAM SIZE = {} mm'.format(self.SIM.beamsize_mm))
    print('SIMTBX: EXPOSURE = {} SEC'.format(self.SIM.exposure_s))

    # now actually run the simulation
    print('SIMTBX: ADDING SPOTS...')
    self.SIM.add_nanoBragg_spots()

    if self.crystal_size is not None:
      scale = (self.crystal_size/1000)**3 / (self.SIM.xtal_size_mm[0] *
                                              self.SIM.xtal_size_mm[1] *
                                              self.SIM.xtal_size_mm[2])
      print('SIMTBX: SCALING BY {} '.format(scale))

      self.SIM.raw_pixels *= scale

    if self.add_background:
      print('SIMTBX: ADDING WATER BACKGROUND')
      # rough approximation to water: interpolation points for sin(theta/lambda) vs structure factor
      bg = flex.vec2_double(
        [(0, 2.57), (0.0365, 2.58), (0.07, 2.8), (0.12, 5), (0.162, 8),
         (0.2, 6.75), (0.18, 7.32), (0.216, 6.75), (0.236, 6.5), (0.28, 4.5),
         (0.3, 4.3), (0.345, 4.36), (0.436, 3.77), (0.5, 3.17)])
      self.SIM.Fbg_vs_stol = bg
      self.SIM.amorphous_sample_thick_mm = 0.1
      self.SIM.amorphous_density_gcm3 = 1
      self.SIM.amorphous_molecular_weight_Da = 18
      self.SIM.add_background()

      print('SIMTBX: ADDING AIR BACKGROUND')
      bg = flex.vec2_double(
        [(0, 14.1), (0.045, 13.5), (0.174, 8.35), (0.35, 4.78), (0.5, 4.22)])
      self.SIM.Fbg_vs_stol = bg
      self.SIM.amorphous_sample_thick_mm = 35  # between beamstop and collimator
      self.SIM.amorphous_density_gcm3 = 1.2e-3
      self.SIM.amorphous_sample_molecular_weight_Da = 28  # nitrogen = N2
      self.SIM.add_background()

    # set this to 0 or -1 to trigger automatic radius.  could be very slow with bright images
    self.SIM.detector_psf_kernel_radius_pixels = 5
    self.SIM.detector_psf_fwhm_mm = 0.08
    self.SIM.detector_psf_type = shapetype.Fiber

    if self.add_noise:
      print('SIMTBX: ADDING NOISE...')
      self.SIM.add_noise()

    # write out a file on arbitrary scale, header contains beam center in various conventions
    print('DEBUG: FILENAME = ', self.img_file)
    print('SIMTBX: CONVERTING TO SMV FORMAT...')
    self.SIM.to_smv_format(fileout="intimage_001.img")

    # output image file
    if self.img_file is not None:
      self.SIM.to_smv_format(fileout=self.img_file)
      print('IMAGE WRITTEN TO ', self.img_file)

    # send pixels back to GUI when done
    try:
      evt = nanoBraggAllDone(tp_EVT_NBDONE, -1, pixels=self.SIM.raw_pixels)
      wx.PostEvent(self.parent, evt)
    except Exception as e:
      print(e)

  def fcalc_from_pdb(self, resolution=None, algorithm=None, wavelength=0.9):
    ''' Generate FCalc from PDB-formatted coordinates '''

    with open(self.params.reference_coordinates, 'r') as pdb_file:
      pdb_lines = pdb_file.readlines()

    # Read in coordinates
    pdb_inp = pdb.input(source_info=None, lines=pdb_lines)
    xray_structure = pdb_inp.xray_structure_simple()

    # take a detour to insist on calculating anomalous contribution of every atom
    scatterers = xray_structure.scatterers()
    for sc in scatterers:
      expected_henke = henke.table(sc.element_symbol()).at_angstrom(wavelength)
      sc.fp = expected_henke.fp()
      sc.fdp = expected_henke.fdp()

    # how do we do bulk solvent?
    primitive_xray_structure = xray_structure.primitive_setting()
    P1_primitive_xray_structure = primitive_xray_structure.expand_to_p1()
    fcalc = P1_primitive_xray_structure.structure_factors(
      d_min=resolution, anomalous_flag=True, algorithm=algorithm).f_calc()
    return fcalc.amplitudes()

# ------------------------------- MISCELLANEOUS ------------------------------ #

class ImageViewerThread(Thread):
  ''' Worker thread that will move the image viewer launch away from the GUI
  and hopefully will prevent the image selection dialog freezing on MacOS'''
  def __init__(self,
               parent,
               file_string,
               viewer='dials.image_viewer',
               img_type=None):
    Thread.__init__(self)
    self.parent = parent
    self.file_string = file_string
    self.viewer = viewer
    self.img_type = img_type

  def run(self):
    command = '{} {}'.format(self.viewer, self.file_string)
    easy_run.fully_buffered(command)
