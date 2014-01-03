# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# $Id:

"""Base class for dark subtraction, gain masking and creation of dxtbx ojbect
"""
from __future__ import division

__version__ = ""

from scitbx.array_family import flex
from xfel.cxi.cspad_ana import cspad_tbx
from xfel.cxi.cspad_ana.mod_event_info import mod_event_info
from xfel.cftbx.detector import cspad_cbf_tbx
import os, copy
from dxtbx.format.Registry import Registry
import numpy as np

class mod_cspad_cbf(mod_event_info):
  """Dark subtraction and gain masking
  """


  def __init__(self,
               address,
               metrology=None,
               dark_path=None,
               mask_path=None,
               gain_map_path=None,
               gain_map_level=None,
               **kwds):
    """The mod_cspad_cbf class constructor stores the
    parameters passed from the pyana configuration file in instance
    variables.

    @param address         Address string XXX Que?!
    @param metrology       CBF header containing metrology info
    @param dark_path       Path to input average dark image
    @param mask_path       Path to input mask.  Pixels to mask out should be set to -2
    @param gain_map_path   Path to input gain map.  Multplied times the image.
    @param gain_map_level  If set, all the '1' pixels in the gain_map are set to this multiplier
                           and all the '0' pixels in the gain_map are set to '1'. If not set,
                           use the values in the gain_map directly
    """

    super(mod_cspad_cbf, self).__init__(address=address, **kwds)

    self.dark_path = cspad_tbx.getOptEvalOrString(dark_path)
    self.mask_path = cspad_tbx.getOptEvalOrString(mask_path)
    gain_map_path = cspad_tbx.getOptString(gain_map_path)
    self.gain_map_level = cspad_tbx.getOptFloat(gain_map_level)

    self.cspad_img = None # The current image - set by self.event()

    # Get and parse metrology.
    device = cspad_tbx.address_split(self.address)[2]
    assert device == 'Cspad'
    self.metrology = cspad_tbx.getOptString(metrology)
    assert os.path.isfile(self.metrology)
    self.reader = Registry.find(self.metrology)


    # Load the dark image and ensure it is signed and at least 32 bits
    # wide, since it will be used for differencing.  If a dark image
    # is provided, a standard deviation image is required, and all the
    # ADU scales must match up.
    self.dark_img = None
    if (self.dark_path is not None):
      reader = Registry.find(self.dark_path)
      self.dark_img = reader(self.dark_path)

    # Load the mask image and ensure it is signed and at least 32 bits
    # wide, since it will be used for differencing.
    self.mask_img = None
    if (self.mask_path is not None):
      pass
      #mask_dict = easy_pickle.load(self.mask_path)
      #self.mask_img = mask_dict['DATA']
      #assert isinstance(self.mask_img, flex.double) or isinstance(self.mask_img, flex.int)

    self.gain_map = None
    if gain_map_path is not None:
      reader = Registry.find(gain_map_path)
      gain_img = reader(gain_map_path)
      self.gain_map = [gain_img.get_raw_data(i) for i in xrange(len(gain_img.get_detector()))]
      if self.gain_map_level is not None:
        cspad_tbx.dynamic_range *= self.gain_map_level
        for i, map in enumerate(self.gain_map):
          sel = flex.bool([bool(f) for f in map])
          sel.reshape(flex.grid(map.focus()))
          map = map.set_selected(~sel, self.gain_map_level)
          self.gain_map[i] = map.set_selected(sel, 1)
      assert False not in [isinstance(map, flex.double) for map in self.gain_map]


  def beginjob(self, evt, env):
    """The beginjob() function does one-time initialisation from
    event- or environment data.  It is called at an XTC configure
    transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """

    reader = Registry.find(self.metrology)
    self.base_dxtbx = reader(self.metrology)
    self.base_cbf = self.base_dxtbx._cbf_handle
    self.base_dxtbx._cbf_handle = None # have to do this because copy.deep_copy fails on
                                       # swig objects, which _cbf_handle is based on

    super(mod_cspad_cbf, self).beginjob(evt, env)

  def event(self, evt, env):
    """The event() function is called for every L1Accept transition.

    @param evt Event data object, a configure object
    @param env Environment object
    """
    super(mod_cspad_cbf, self).event(evt, env)
    if (evt.get("skip_event")):
      return

    # This module only applies to detectors for which a distance is
    # available.
    distance = cspad_tbx.env_distance(self.address, env, self._detz_offset)
    if distance is None:
      self.nfail += 1
      self.logger.warning("event(): no distance, shot skipped")
      evt.put(True, "skip_event")
      return

    # Early return if the full detector image is already stored in the
    # event.  Otherwise, get it from the stream as a double-precision
    # floating-point flex array.
    self.cspad_img = evt.get(self.address)
    if self.cspad_img is None:
      # Set up the dxtbx object and its cbf handle
      self.cspad_img = copy.deepcopy(self.base_dxtbx)
      self.cspad_img._cbf_handle = cbf = cspad_cbf_tbx.copy_cbf_header(self.base_cbf)
      cbf.set_datablockname("cspad_" + self.timestamp)

      quads = evt.getCsPadQuads(self.address, env)
      n_asics = sum([len(quad.data()) for quad in quads]) * 2
      cspad_cbf_tbx.add_frame_specific_cbf_tables(cbf, self.wavelength,self.timestamp,
        [cspad_tbx.dynamic_range]*n_asics)

      # Set the distance, I.E., the length translated along the Z axis
      cbf.find_category("diffrn_scan_frame_axis")
      cbf.find_column("axis_id")
      cbf.find_row("AXIS_D0_Z") # XXX discover the Z axis somehow, don't use D0 here
      cbf.find_column("displacement")
      cbf.set_value(str(-distance))

      # Explicitly reset the detector object now that the distance is set correctly
      self.cspad_img._detector_instance = self.cspad_img._detector()

      # Explicitly set up the beam object now that the tables are all loaded correctly
      self.cspad_img._beam_instance = self.cspad_img._beam()

      # Get the data and add it to the cbf handle
      tiles = {}
      for quad in quads:
        for s in xrange(len(quad.data())):
          sensor = quad.data()[s]
          for a in xrange(2):
            tiles[(0,quad.quad(),s,a)] = flex.double(sensor[:,194*a:194*(a+1)].astype(np.float64))

      # If a dark image was provided, subtract it from the image.
      if (self.dark_img is not None):
        assert len(tiles) == len(self.dark_img.get_detector())
        for i, k in enumerate(sorted(tiles)):
          tiles[k] -= self.dark_img.get_raw_data(i)

      # If a gain map was provided, multiply it times the image
      if self.gain_map is not None:
        assert len(tiles) == len(self.gain_map)
        for i, k in enumerate(sorted(tiles)):
          tiles[k] *= self.gain_map[i]

      # add the pixel data
      cspad_cbf_tbx.add_tiles_to_cbf(cbf,tiles)

      self.logger.info("Processed %s"%self.timestamp)
    else:
      #if (self.mask_img is not None):
        #sel = self.mask_img == -2
        #self.cspad_img.set_selected(sel, -2)
      return

    if self.cspad_img is None:
      self.nfail += 1
      self.logger.warning("event(): no image, shot skipped")
      evt.put(True, "skip_event")
      return

    #if (self.mask_img is not None):
    #  sel = self.mask_img == -2
    #  self.cspad_img.set_selected(sel, -2)


  def endjob(self, env):
    """
    @param env Environment object
    """

    if 0 and self.dark_path is not None and self.nmemb > 1:
      pass
